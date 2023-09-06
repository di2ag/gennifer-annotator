### Gunicorn Configuration File ###

timeout = 360
graceful_timeout = 360
limit_request_field_size  = 0
limit_request_line = 0
limit_request_fields = 0
proxy_allow_ips = '*'
workers=2
errorlog='gunicorn-error.log'
accesslog='gunicorn-access.log'
loglevel='debug'
