# this will take a very long time! don't actually do it!
for file in network_setup/transformations_set_1/*.json; do
  relpath=${file:36}  # strip off "network_setup/transformations/"
  dirpath=${relpath%.*}  # strip off final ".json"
  # loop over three repeats
  for repeat in {1..3}; do
      openfe quickrun $file -o results/repeat${repeat}/$relpath -d results/repeat${repeat}/$dirpath
  done
done