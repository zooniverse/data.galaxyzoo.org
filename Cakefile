{ spawn } = require 'child_process'

run = (fullCommand) ->
  [command, args...] = fullCommand.split /\s+/
  child = spawn command, args
  child.stdout.on 'data', process.stdout.write.bind process.stdout
  child.stderr.on 'data', process.stderr.write.bind process.stderr

task 'watch-stylus', 'Recompile Stylus files when they change', ->
  console.log 'Watching .styl files in ./styles'
  run 'stylus -u ./node_modules/nib/lib/nib -w ./styles -o ./public/styles'

task 'watch-coffee', 'Watch CoffeeScript changes during development', ->
  console.log 'Watching .coffee files in ./scripts'
  run 'coffee -o ./public/scripts -cw ./scripts'

task 'compile', 'Compile both Stylus and Coffeescript', ->
  run 'stylus -u ./node_modules/nib/lib/nib -c ./styles -o ./public/styles'
  run 'coffee -o ./public/scripts -c ./scripts'

task 'serve', 'Run a dev server', ->
  invoke 'watch-coffee' 
  invoke 'watch-stylus'

  console.log "Running a server at localhost:8080"
  run 'static ./public'
