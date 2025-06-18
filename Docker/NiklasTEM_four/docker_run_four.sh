# same command as above with /home/nab/Niklas/DockerShare/DockerTEMNiklasThree as a base path

sudo docker stop neo4j-niklas-tem-niklas-four

sudo docker remove neo4j-niklas-tem-niklas-four

sudo docker run -d --name neo4j-niklas-tem-niklas-four \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/niklasniklaspwtemnewstart \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_PLUGINS='["apoc", "graph-data-science"]' \
  -e NEO4J_dbms_security_procedures_unrestricted="apoc.*,gds.*" \
  -e NEO4J_dbms_security_procedures_allowlist="apoc.*,gds.*" \
  -e NEO4JLABS_PLUGINS='["apoc", "graph-data-science"]' \
  -p 2147:7687 \
  -p 2148:7474 \
  -v /home/nab/Niklas/DockerShare/DockerTEMNiklasFour/data:/data \
  -v /home/nab/Niklas/DockerShare/DockerTEMNiklasFour/import:/import \
  -v /home/nab/Niklas/DockerShare/DockerTEMNiklasFour/plugins:/plugins \
  my-neo4j:latest