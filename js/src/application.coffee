
TopBar = window.zooniverse.controllers.TopBar
(new TopBar).el.prependTo 'body'

Footer = window.zooniverse.controllers.Footer
(new Footer).el.appendTo 'body > .footer > .container'