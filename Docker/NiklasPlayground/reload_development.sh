sudo docker stop pyeed-neo4j-niklas
sudo docker remove pyeed-neo4j-niklas
sudo docker build --no-cache -t pyeed-neo4j-niklas .
sudo docker run -d --name pyeed-neo4j-niklas \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/12345678 \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_dbms_unmanaged__extension__classes="n10s.endpoint=/rdf" \
  -e NEO4J_dbms_security_procedures_unrestricted='["n10s.*", "apoc.*"]' \
  -e NEO4J_PLUGINS='["graph-data-science", "apoc"]' \
  -e NEO4J_dbms_security_procedures_allowlist='["n10s.*", "apoc.*"]' \
  -p 7474:7474 \
  -p 7687:7687 \
  -v $(pwd)/data:/data \
  -v $(pwd)/import:/import \
  neo4j:latest