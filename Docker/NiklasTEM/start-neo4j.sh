#!/bin/bash

# Ensure any failure in the script stops the script
echo "Starting Neo4j with NEO4J_AUTH set to: $NEO4J_AUTH"
# Start Neo4j in the foreground
exec neo4j console