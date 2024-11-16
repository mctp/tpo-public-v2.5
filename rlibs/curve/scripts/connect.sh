#!/bin/bash

if [[ -n "$CURVE_DB" ]]; then
  args=$(echo $CURVE_DB | sed 's_:_ _g' | xargs printf "host=%s port=%s dbname=%s user=%s password=%s")
  echo "psql $args"
  psql "$args"
else
  echo "please set CURVE_DB env variable"
fi
