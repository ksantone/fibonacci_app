FROM python:3
ENV PYTHONUNBUFFERED=1
COPY . .
RUN pip install -r requirements.txt
RUN cd pipelines/MyBindings/build
RUN apt-get update -y
RUN apt-get update
RUN apt-get install cmake
RUN rm CMakeCache.txt
RUN cmake ..
RUN make

CMD [ "python", "./manage.py", "runserver", "0.0.0.0:8000" ]
