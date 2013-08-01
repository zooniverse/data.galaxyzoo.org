{spawn} = require 'child_process'

run = (fullCommand) ->
  [command, args...] = fullCommand.split /\s+/
  child = spawn command, args
  child.stdout.on 'data', process.stdout.write.bind process.stdout
  child.stderr.on 'data', process.stderr.write.bind process.stderr

task 'watch-stylus', 'Recompile Stylus files when they change', ->
  console.log 'Watching .styl files in ./stylesheets/src'
  run 'stylus -u ./node_modules/nib/lib/nib -w ./stylesheets/src -o ./stylesheets'

task 'watch-coffee', 'Watch CoffeeScript changes during development', ->
  console.log 'Watching for CoffeeScript in ./src'
  run 'coffee -w -o ./js -c ./js/src'

task 'watch-jade', 'Watch jade template files', ->
  console.log 'Watching .jade files in ./templates'
  run 'jade -P -o ./ -w ./templates/index.jade' # lol

task 'serve', 'Run a dev server', ->
  invoke 'watch-coffee'
  invoke 'watch-stylus'
  invoke 'watch-jade'

  console.log "Running a server at localhost:8080"
  run 'static'

