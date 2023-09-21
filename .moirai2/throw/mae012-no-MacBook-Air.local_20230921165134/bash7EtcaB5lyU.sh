#!/bin/bash
cd /Users/ah3q/Sites/github/HDRGenome
if [[ "$(/usr/local/bin/docker images -q moirai2/biotools 2> /dev/null)" == "" ]]; then
echo "'moirai2/biotools' docker doesn't exist" > .moirai2/202309211651340000_f861d2f9b1542f890f9d3238bb6dfba640526_local_1/stderr.txt
echo "error	"`date +%s` > .moirai2/202309211651340000_f861d2f9b1542f890f9d3238bb6dfba640526_local_1/status.txt
touch .moirai2/202309211651340000_f861d2f9b1542f890f9d3238bb6dfba640526_local_1/stdout.txt
touch .moirai2/202309211651340000_f861d2f9b1542f890f9d3238bb6dfba640526_local_1/log.txt
exit
else
/usr/local/bin/docker \
  run \
  --rm \
  --workdir=/root \
  -u `id -u`:`id -g` \
  -v '/Users/ah3q/Sites/github/HDRGenome:/root' \
  moirai2/biotools \
  /bin/bash /root/.moirai2/202309211651340000_f861d2f9b1542f890f9d3238bb6dfba640526_local_1/run.sh \
  > .moirai2/202309211651340000_f861d2f9b1542f890f9d3238bb6dfba640526_local_1/stdout.txt \
  2> .moirai2/202309211651340000_f861d2f9b1542f890f9d3238bb6dfba640526_local_1/stderr.txt
fi
if [ -e .moirai2/throw/mae012-no-MacBook-Air.local_20230921165134/bash7EtcaB5lyU.stdout ] && [ ! -s .moirai2/throw/mae012-no-MacBook-Air.local_20230921165134/bash7EtcaB5lyU.stdout ]; then
rm .moirai2/throw/mae012-no-MacBook-Air.local_20230921165134/bash7EtcaB5lyU.stdout
fi
if [ -e .moirai2/throw/mae012-no-MacBook-Air.local_20230921165134/bash7EtcaB5lyU.stderr ] && [ ! -s .moirai2/throw/mae012-no-MacBook-Air.local_20230921165134/bash7EtcaB5lyU.stderr ]; then
rm .moirai2/throw/mae012-no-MacBook-Air.local_20230921165134/bash7EtcaB5lyU.stderr
fi
rm /Users/ah3q/Sites/github/HDRGenome/.moirai2/throw/mae012-no-MacBook-Air.local_20230921165134/bash7EtcaB5lyU.sh
