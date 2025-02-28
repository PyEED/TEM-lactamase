# same command as above with /home/nab/Niklas/DockerShare/DockerTEMClean as a base path

sudo docker stop neo4j-niklas-tem-clean

sudo docker remove neo4j-niklas-tem-clean

sudo docker run -d --name neo4j-niklas-tem-clean \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/niklasniklaspwtemclean \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_PLUGINS='["apoc", "graph-data-science"]' \
  -e NEO4J_dbms_security_procedures_unrestricted="apoc.*,gds.*" \
  -e NEO4J_dbms_security_procedures_allowlist="apoc.*,gds.*" \
  -e NEO4JLABS_PLUGINS='["apoc", "graph-data-science"]' \
  -p 2123:7687 \
  -p 2124:7474 \
  -v /home/nab/Niklas/DockerShare/DockerTEMClean/data:/data \
  -v /home/nab/Niklas/DockerShare/DockerTEMClean/import:/import \
  -v /home/nab/Niklas/DockerShare/DockerTEMClean/plugins:/plugins \
  my-neo4j:latest
