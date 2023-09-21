#!/bin/bash
cd /Users/ah3q/Sites/github/HDRGenome
if [[ "$(/usr/local/bin/docker images -q moirai2/biotools 2> /dev/null)" == "" ]]; then
echo "'moirai2/biotools' docker doesn't exist" > .moirai2/202309211651270000_0d01f418161c755231926823aeac1907679_local_1/stderr.txt
echo "error	"`date +%s` > .moirai2/202309211651270000_0d01f418161c755231926823aeac1907679_local_1/status.txt
touch .moirai2/202309211651270000_0d01f418161c755231926823aeac1907679_local_1/stdout.txt
touch .moirai2/202309211651270000_0d01f418161c755231926823aeac1907679_local_1/log.txt
exit
else
/usr/local/bin/docker \
  run \
  --rm \
  --workdir=/root \
  -u `id -u`:`id -g` \
  -v '/Users/ah3q/Sites/github/HDRGenome:/root' \
  moirai2/biotools \
  /bin/bash /root/.moirai2/202309211651270000_0d01f418161c755231926823aeac1907679_local_1/run.sh \
  > .moirai2/202309211651270000_0d01f418161c755231926823aeac1907679_local_1/stdout.txt \
  2> .moirai2/202309211651270000_0d01f418161c755231926823aeac1907679_local_1/stderr.txt
fi
if [ -e .moirai2/throw/mae012-no-MacBook-Air.local_20230921165127/bashxikopEKMOZ.stdout ] && [ ! -s .moirai2/throw/mae012-no-MacBook-Air.local_20230921165127/bashxikopEKMOZ.stdout ]; then
rm .moirai2/throw/mae012-no-MacBook-Air.local_20230921165127/bashxikopEKMOZ.stdout
fi
if [ -e .moirai2/throw/mae012-no-MacBook-Air.local_20230921165127/bashxikopEKMOZ.stderr ] && [ ! -s .moirai2/throw/mae012-no-MacBook-Air.local_20230921165127/bashxikopEKMOZ.stderr ]; then
rm .moirai2/throw/mae012-no-MacBook-Air.local_20230921165127/bashxikopEKMOZ.stderr
fi
rm /Users/ah3q/Sites/github/HDRGenome/.moirai2/throw/mae012-no-MacBook-Air.local_20230921165127/bashxikopEKMOZ.sh
