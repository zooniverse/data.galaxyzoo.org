data.galaxyzoo.org
====================

To run locally, first make sure the repo files are updated. Then run the command:

    python -m SimpleHTTPServer 8000

Then open the following in the browser to view the page:

    localhost:8000/public/

Other static file servers can also be run in the root directory if you prefer. 

To deploy changes, use node v0.10.x, eg:

```
nvm use v0.10
npm install
npm run-script deploy
```

Make sure you have variables set for `AMAZON_ACCESS_KEY_ID` and `AMAZON_SECRET_ACCESS_KEY` before running the final deploy step.
