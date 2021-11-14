FROM python:3
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y cmake

WORkDIR /app

COPY . .

RUN pip install -r requirements.txt

RUN mkdir pipelines/BioLib/build \
    && cd pipelines/BioLib/build \
    && cmake .. \
    && make \
    && mv denovo_sequencing.*.so /app/denovo_sequencing.so

ENV PYTHONPATH=/app

CMD [ "python", "./manage.py", "runserver", "0.0.0.0:8000" ]
