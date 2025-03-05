import os
import sys
from obspy.clients.filesystem.tsindex import Indexer

# for this example get the file path to test data
filepath = sys.argv[1]
# filepath = 'test/'
# filepath = 'waveforms/'
database = os.path.join(filepath, 'timeseries.sqlite')

# create a new Indexer instance
indexer = Indexer(filepath, database=database, leap_seconds_file='leap-seconds.list', filename_pattern='*.mseed')

indexer.run(build_summary=True, relative_paths=True, reindex=True)
