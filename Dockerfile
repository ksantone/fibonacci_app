FROM python:3
ENV PYTHONUNBUFFERED=1
COPY . .
RUN pip install -r requirements.txt
RUN cd pipelines/MyBindings/build
RUN curl https://cmake.org/files/v3.12/cmake-3.12.0-Linux-x86_64.sh -o cmake.sh
RUN sh cmake.sh --prefix=${pwd}/cmake-test-install
RUN cmake ..
RUN make

CMD [ "python", "./manage.py", "runserver", "0.0.0.0:8000" ]
