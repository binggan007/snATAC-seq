import request
import sys

exp = sys.argv[1]
file = sys.argv[2]
usr = sys.argv[3]
pswd = sys.argv[4]

auth = (usr,pswd)
r = requests.get('https://www.encodeproject.org/files/'+exp+'/@@download/'+exp+'.tar.gz', auth=auth)
open(exp+'.tar.gz', 'wb').write(r.content)
