data.galaxyzoo.org
====================

To run locally, first make sure the repo files are updated. Then run the command:

    python -m SimpleHTTPServer 8000
    
for Python 2.*, or

    python -m http.server 8000
    
for Python 3.*.

Then open the following in the browser to view the page:

    localhost:8000/public/

Other static file servers can also be run in the root directory if you prefer. 

Deploy changes with Docker:

```
docker-compose build
docker-compose run dev npm run-script deploy
```

Alternatively, use node v9, eg:

```
nvm use v9
npm install
npm run-script deploy
```

Make sure you have variables set for `AMAZON_ACCESS_KEY_ID` and `AMAZON_SECRET_ACCESS_KEY` before running the final deploy step.
