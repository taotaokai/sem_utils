from obspy.clients.filesystem.tsindex import Client
import os

# for this example get the file path to test data
filepath = get_test_data_filepath()


db_name = 'timeseries.sqlite'
data_dir = 'waveforms/'
station_dir = 'stations/'

# create a new Client instance
client = Client(db_name, datapath_replace=("^", waveforms))

client
