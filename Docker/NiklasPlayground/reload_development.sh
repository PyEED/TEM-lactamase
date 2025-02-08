sudo docker stop neo4j-niklas-playground
sudo docker remove neo4j-niklas-playground
sudo docker build --no-cache -t neo4j-niklas-playground .
sudo docker run -d --name neo4j-niklas-playground \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/12345678 \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_dbms_unmanaged__extension__classes="n10s.endpoint=/rdf" \
  -e NEO4J_dbms_security_procedures_unrestricted='["n10s.*", "apoc.*"]' \
  -e NEO4J_PLUGINS='["graph-data-science", "apoc"]' \
  -e NEO4J_dbms_security_procedures_allowlist='["n10s.*", "apoc.*"]' \
  -p 7474:7474 \
  -p 7687:7687 \
  -v /home/nab/Niklas/DockerShare/DockerPlayground/data:/data \
  -v /home/nab/Niklas/DockerShare/DockerPlayground/import:/import \
  -v /home/nab/Niklas/DockerShare/DockerPlayground/plugins:/plugins \
  my-neo4j:latest