import copy
from typing import List, Dict, Tuple, Optional, Union
from obspy.core.inventory import Inventory, Network, Station, Channel

def isolate_and_merge_station(
    inv: Inventory, 
    network_id: str, 
    station_id: str, 
    ignore_station_creation_date: bool = False
) -> Inventory:
    """
    Takes an inventory object, isolates the given station and merges them.

    Merging is sometimes necessary as many files have the same station
    multiple times.

    Returns the processed inventory object. The original one will not be
    changed.

    :param inv: The inventory.
    :param network_id: The network id.
    :param station_id: The station id.
    :param ignore_station_creation_date: Ignore the creation date of the station when merging.
    :return: Processed inventory object with merged stations and channels.
    """
    # Create a deep copy and isolate the requested network/station
    inv = copy.deepcopy(
        inv.select(network=network_id, station=station_id, keep_empty=True)
    )

    # Handle case where no networks match the selection
    if len(inv.networks) == 0:
        return inv

    # Merge networks if necessary - it should be safe to always merge networks
    # with the same code
    if len(inv.networks) > 1:
        inv = _merge_networks(inv)

    # Merge stations if possible
    if len(inv.networks[0].stations) > 1:
        inv.networks[0].stations = _merge_stations(
            inv.networks[0].stations, 
            ignore_station_creation_date
        )

    # Merge channels if possible
    for station in inv.networks[0].stations:
        # Only attempt to merge if there's more than one channel
        if len(station.channels) > 1:
            station.channels = _merge_channels(station.channels)

    return inv

def _merge_networks(inv: Inventory) -> Inventory:
    """
    Merge multiple networks with the same code into a single network.
    
    :param inv: Inventory containing networks to merge
    :return: Inventory with merged networks
    """
    network = inv.networks[0]
    for other_network in inv.networks[1:]:
        # Merge the stations
        network.stations.extend(other_network.stations)
        
        # Update the times if necessary
        if other_network.start_date is not None:
            if (network.start_date is None or 
                network.start_date > other_network.start_date):
                network.start_date = other_network.start_date
                
        # None is the "biggest" end_date
        if network.end_date is not None and other_network.end_date is not None:
            if other_network.end_date > network.end_date:
                network.end_date = other_network.end_date
        elif other_network.end_date is None:
            network.end_date = None
            
        # Update comments
        network.comments = list(
            set(network.comments).union(set(other_network.comments))
        )
        
        # Update the number of stations
        if other_network.total_number_of_stations:
            if (not network.total_number_of_stations or 
                network.total_number_of_stations < other_network.total_number_of_stations):
                network.total_number_of_stations = other_network.total_number_of_stations
                
        # Update the other elements
        network.alternate_code = network.alternate_code or other_network.alternate_code or None
        network.description = network.description or other_network.description or None
        network.historical_code = network.historical_code or other_network.historical_code or None
        network.restricted_status = network.restricted_status or other_network.restricted_status
        
    inv.networks = [network]
    return inv

def _check_stations_meta_equal(
    station_a: Station, 
    station_b: Station, 
    ignore_station_creation_date: bool
) -> bool:
    """
    Check if two stations have equal metadata (excluding channels, dates, and optionally creation date).
    
    :param station_a: First station to compare
    :param station_b: Second station to compare
    :param ignore_station_creation_date: Whether to ignore creation date in comparison
    :return: True if stations are equal, False otherwise
    """
    ignored_keys = ["channels", "start_date", "end_date"]
    # Whether to ignore the creation date
    if ignore_station_creation_date:
        ignored_keys.append("_creation_date")
        
    # Shallow copies
    content_a = copy.copy(station_a.__dict__)
    content_b = copy.copy(station_b.__dict__)
    
    for content in [content_a, content_b]:
        # Get rid of keys we don't want to test
        for key in ignored_keys:
            if key in content:
                del content[key]
                
    return content_a == content_b

def _merge_stations(stations: List[Station], ignore_station_creation_date: bool) -> List[Station]:
    """
    Merge stations with identical metadata into single stations with extended time ranges.
    
    :param stations: List of stations to merge
    :param ignore_station_creation_date: Whether to ignore creation date when comparing stations
    :return: List of merged stations
    """
    # Always keep the first one
    merged_stations = [stations[0]]
    
    # Now loop over others, see if they match an existing one and add the
    # channels while adjusting start and end times
    for new_station in stations[1:]:
        for existing_station in merged_stations:
            if not _check_stations_meta_equal(
                existing_station, 
                new_station, 
                ignore_station_creation_date
            ):
                continue
                
            # Here we have a station where everything except the channels
            # and start and end date are identical. Merge both.
            _merge_station_time_ranges(existing_station, new_station)
            
            # Merge the channels
            existing_station.channels.extend(new_station.channels)
            
            # Finally break the loop because the station has been dealt with
            break
        else:
            merged_stations.append(new_station)
            
    return merged_stations

def _merge_station_time_ranges(existing_station: Station, new_station: Station) -> None:
    """
    Update time ranges for an existing station based on a new station's time ranges.
    
    :param existing_station: Station to update
    :param new_station: Station with new time ranges
    """
    # Update the start date if the start date from the new station is earlier
    # than the existing start date or if there is no existing start date
    if new_station.start_date is not None:
        if (existing_station.start_date is None or 
            existing_station.start_date > new_station.start_date):
            existing_station.start_date = new_station.start_date
            
    # None is the "biggest" end_date
    if existing_station.end_date is not None and new_station.end_date is not None:
        if new_station.end_date > existing_station.end_date:
            existing_station.end_date = new_station.end_date
    elif new_station.end_date is None:
        existing_station.end_date = None

def _check_channels_meta_equal(channel_a: Channel, channel_b: Channel) -> bool:
    """
    Check if two channels have equal metadata (excluding dates).
    
    :param channel_a: First channel to compare
    :param channel_b: Second channel to compare
    :return: True if channels are equal, False otherwise
    """
    ignored_keys = ["start_date", "end_date"]
    # Shallow copies
    content_a = copy.copy(channel_a.__dict__)
    content_b = copy.copy(channel_b.__dict__)
    
    for content in [content_a, content_b]:
        # Get rid of keys we don't want to test
        for key in ignored_keys:
            if key in content:
                del content[key]
                
    return content_a == content_b

def _merge_channels(channels: List[Channel]) -> List[Channel]:
    """
    Merge channels with identical metadata into single channels with extended time ranges.
    
    :param channels: List of channels to merge
    :return: List of merged channels
    """
    if not channels:
        return channels
        
    # Group channels by their metadata
    channel_groups: List[Dict] = [{
        "channel": channels[0],
        "time_ranges": [(channels[0].start_date, channels[0].end_date)],
    }]
    
    for channel in channels[1:]:
        for group in channel_groups:
            if _check_channels_meta_equal(group["channel"], channel):
                group["time_ranges"].append((channel.start_date, channel.end_date))
                break
        else:
            new_group = {
                "channel": channel,
                "time_ranges": [(channel.start_date, channel.end_date)],
            }
            channel_groups.append(new_group)

    # Create unique channels with merged time ranges
    unique_channels = []
    for group in channel_groups:
        channel = group["channel"]
        start_dates = [t[0] for t in group["time_ranges"] if t[0]]
        end_dates = [t[1] for t in group["time_ranges"] if t[1]]
        
        channel.start_date = min(start_dates) if start_dates else None
        channel.end_date = max(end_dates) if end_dates else None
        unique_channels.append(channel)

    return unique_channels