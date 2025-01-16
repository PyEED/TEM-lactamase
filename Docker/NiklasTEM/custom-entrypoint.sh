#!/bin/bash
# Set up a trap to catch when the Neo4j service stops and take necessary action
trap "echo 'Neo4j has stopped. Container will stay up.'; sleep infinity" SIGTERM SIGINT

# Run Neo4j in the foreground
neo4j console &
NEO4J_PID=$!

# Wait for Neo4j process to exit
wait $NEO4J_PID

# In case the Neo4j process exits, keep the container running
sleep infinity
