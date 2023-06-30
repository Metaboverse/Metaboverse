# Use an official Python runtime as a parent image
FROM python:3.9

# Set the working directory in the container to /app
WORKDIR /app

COPY __version__.txt __version__.txt
COPY LICENSE LICENSE
COPY README.md README.md

COPY requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install your application 
COPY . .
RUN pip install -e .

# Make port 80 available to the world outside this container
EXPOSE 80

# Run the application
CMD ["python", "-m", "metaboverse_cli"]
