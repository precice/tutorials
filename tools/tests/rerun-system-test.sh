#!/usr/bin/env sh
set -e -u

cd "$(dirname "$0")"

# Unzipped CI artifacts are often owned by the host user. The prepare container
# runs as precice and must edit precice-config.xml in this directory.
chmod -R a+rwX .

echo "[systemtests] Building tutorial images..."
docker compose --file docker-compose.tutorial.yaml build

echo "[systemtests] Running tutorial containers..."
docker compose --file docker-compose.tutorial.yaml up

if [ -f docker-compose.field_compare.yaml ]; then
  echo "[systemtests] Running fieldcompare..."
  docker compose --file docker-compose.field_compare.yaml up --exit-code-from field-compare
else
  echo "[systemtests] Skipping fieldcompare (docker-compose.field_compare.yaml not present; the original run likely failed before compare)."
fi
