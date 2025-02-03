import os
import subprocess
from datetime import datetime

from dotenv import load_dotenv

# Neo4j Connection Info
URI = "http://129.69.129.130:8124"
USER = "neo4j"
CONTAINER_NAME = "neo4j-niklas-tem"

load_dotenv()
PASSWORD = os.getenv("NEO4J_NIKLAS_TEM_HARRY")
if PASSWORD is None:
    raise ValueError("KEY is not set in the .env file.")


# Backup directory
BACKUP_DIR = "/home/nab/Niklas/DockerShare/DockerTem/backups"
IMPORT_DIR = "/home/nab/Niklas/DockerShare/DockerTem/import"


def run_command(command):
    """Run a shell command and handle errors."""
    try:
        result = subprocess.run(
            command, shell=True, check=True, capture_output=True, text=True
        )
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")
        raise


def backup_neo4j():
    """
    Create a Neo4j database backup.
    """

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    backup_file = f"{BACKUP_DIR}/neo4j_backup_{timestamp}.dump"

    # Stop the database
    print("Stopping Neo4j database...")
    run_command(f"sudo docker exec {CONTAINER_NAME} neo4j-admin server stop")

    # Create backup
    command = f"sudo docker exec {CONTAINER_NAME} neo4j-admin database dump --verbose --to-path=/import/ --overwrite-destination=true neo4j"
    print("Creating backup...")
    run_command(command)

    # Start the database again
    print("Starting Neo4j database...")
    run_command(f"sudo docker exec {CONTAINER_NAME} neo4j-admin server start")

    # Move backup from container's import folder to the host
    host_backup_path = os.path.join(IMPORT_DIR, "neo4j.dump")
    run_command(f"mv {host_backup_path} {backup_file}")

    print(f"Backup saved to: {backup_file}")


def restore_neo4j(backup_file):
    """
    Restore a Neo4j database from a backup.
    """

    if not os.path.exists(backup_file):
        raise FileNotFoundError(f"Backup file not found: {backup_file}")

    print("Stopping Neo4j container...")
    run_command(f"sudo docker stop {CONTAINER_NAME}")

    print("Copying backup file to import directory...")
    run_command(f"cp {backup_file} {IMPORT_DIR}/neo4j.dump")

    print("starting Neo4j container...")
    run_command(f"sudo docker start {CONTAINER_NAME}")

    print("Restoring database...")
    command = f"sudo docker exec {CONTAINER_NAME} neo4j-admin database load --from-path=/import/ --verbose --overwrite-destination=true neo4j"
    run_command(command)

    print("Restarting Neo4j container...")

    run_command(f"sudo docker start {CONTAINER_NAME}")
    print("Restore complete.")


if __name__ == "__main__":
    backup_neo4j()

    # Restore backup file is /mnt/nab/backups/TEM_Development/neo4j_backup_20241130_100940.dump
    # restore_neo4j("/mnt/nab/backups/TEM_Development/neo4j_backup_20241219_133125.dump")
