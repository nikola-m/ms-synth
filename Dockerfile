# Containerised reproduction environment for ms-synth v1.0.0.
# Build: docker build -t ms-synth:1.0.0 .
# Run:   docker run --rm ms-synth:1.0.0            (verify published results)
#        docker run --rm -v "$PWD/out:/app/results" ms-synth:1.0.0 make all
FROM python:3.12-slim
WORKDIR /app
COPY requirements-lock.txt pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir -r requirements-lock.txt && pip install --no-cache-dir --no-deps -e .
COPY . .
CMD ["python", "scripts/verify_reproduction.py"]
