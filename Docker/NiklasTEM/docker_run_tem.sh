
sudo docker stop neo4j-niklas-tem

sudo docker remove neo4j-niklas-tem

sudo docker run -d --name neo4j-niklas-tem \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/niklasniklaspwtem \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_PLUGINS='["apoc", "graph-data-science"]' \
  -e NEO4J_dbms_security_procedures_unrestricted="apoc.*,gds.*" \
  -e NEO4J_dbms_security_procedures_allowlist="apoc.*,gds.*" \
  -e NEO4JLABS_PLUGINS='["apoc", "graph-data-science"]' \
  -p 8123:7687 \
  -p 8124:7474 \
  -v /mnt/nab/NiklasTEM/data:/data \
  -v /mnt/nab/NiklasTEM/import:/import \
  -v /mnt/nab/NiklasTEM/plugins:/plugins \
  my-neo4j:latest





