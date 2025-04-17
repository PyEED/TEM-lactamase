# same command as above with /home/nab/Niklas/DockerShare/DockerTEMClean as a base path

sudo docker stop neo4j-niklas-tem-new-start

sudo docker remove neo4j-niklas-tem-new-start

sudo docker run -d --name neo4j-niklas-tem-new-start \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/niklasniklaspwtemnewstart \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_PLUGINS='["apoc", "graph-data-science"]' \
  -e NEO4J_dbms_security_procedures_unrestricted="apoc.*,gds.*" \
  -e NEO4J_dbms_security_procedures_allowlist="apoc.*,gds.*" \
  -e NEO4JLABS_PLUGINS='["apoc", "graph-data-science"]' \
  -p 2127:7687 \
  -p 2128:7474 \
  -v /home/nab/Niklas/DockerShare/DockerNewStart/data:/data \
  -v /home/nab/Niklas/DockerShare/DockerNewStart/import:/import \
  -v /home/nab/Niklas/DockerShare/DockerNewStart/plugins:/plugins \
  my-neo4j:latest
