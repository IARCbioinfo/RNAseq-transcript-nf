#!/bin/bash
cd ~/project/
git config --global user.email "alcalan@fellows.iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git status
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://hub.docker.com/api/build/v1/source/4c39489f-d75e-43d6-9435-f20038dfb4aa/trigger/cb94b4b1-5519-4f4f-9649-b9a10c470026/call/
