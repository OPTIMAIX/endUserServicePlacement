# endUserServicePlacement
The code in this repository implements the algorithm for end user service placement presented in Deliverable 2.1. It improves workload distribution in telecommunications networks, enhancing network efficiency, user experience, and profitability. The algorithm focuses on services characterized by latency, jitter, or traffic density,


To get help for microservice options, type in a terminal:

wsl docker run --rm micro-nfplaning -h


To run the microservice image, type in a terminal:

wsl docker run --rm -p 9900:9910 micro-nfplaning -l trace

Note that the choice of 9900 was arbitrary.  9910 is the default port
exposed by the microservice within the container, and you can map it to any
available port on your host machine. The remainder of this tutorial will
assume port 9900 is used.


Once the microservice image is started, you can check the status of the
service by visiting the following URL in a web browser:

http://localhost:9900/api/health

If the microservice is ready to receive requests, you should see:
"status: ok".


To locally test the running microservice, we will use the curl command.
The generic form of the curl command to post an HTTP request to the
microservice is:

curl -v -H Content-Type:application/json -d '<json-data>' http://<hostname>:<port>/<archivename>/<functionname>

Where:
<json-data> is the input to the called function in json format
<hostname> is the hostname of the machine running the container
<port> is the hostname port that the microservice is listening to
<archivename> is the archive name of the CTF file deployed in the microservice container
<functionname> is the name of the deployed MATLAB function being called by this request


For example, to call function <functionname> with one output and one numeric input:
curl -v -H Content-Type:application/json -d '{"nargout":1,"rhs":[1]}' http://localhost:9900/NFPSarchive/<functionname>

To call function <functionname> with one output and one character array input:
curl -v -H Content-Type:application/json -d '{"nargout":1,"rhs":["abc"]}' http://localhost:9900/NFPSarchive/<functionname>

On Windows, the syntax for the curl command is different.
To call function <functionname> on Windows with one output and one numeric input:
curl -v -H Content-Type:application/json -d "{\"nargout\":1,\"rhs\":[1]}" http://localhost:9900/NFPSarchive/<functionname>

To call function <functionname> on Windows with one output and one character array input:
curl -v -H Content-Type:application/json -d "{\"nargout\":1,\"rhs\":[\"abc\"]}" http://localhost:9900/NFPSarchive/<functionname>

