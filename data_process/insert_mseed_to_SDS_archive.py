import sys
# import sds
from mass_downloader_SDS import sds

mseed_list = sys.argv[1]
archive_root = sys.argv[2]

with open(mseed_list, 'r') as f:
    mseed_files = [l.split()[0] for l in f.readlines() if not l.startswith('#')]

sds_client = sds.Client(archive_root)

for mseed_file in mseed_files:
    print(f'================ {mseed_file}')
    sds_client.insert_mseed_file(mseed_file)
