sudo docker stop neo4j-niklas-tem-main

sudo docker remove neo4j-niklas-tem-main






sudo docker run -d --name neo4j-niklas-tem-main \
  -p 2123:7687 \
  -p 2124:7474 \
  -e NEO4J_AUTH=neo4j/1234567890 \
  --user $(id -u):$(id -g) \
  -v $(pwd)/data:/data \
  -v $(pwd)/import:/import \
  -v $(pwd)/plugins:/plugins \
  my-neo4j:latest
