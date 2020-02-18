(function ($) {
  $.extend($, {
    cupEvent: {
      version: '1.0.0',
      initialized: false,
      uuid: new Date().getTime()
    }
  });

  function Events() {
    this.subscribers = []; //To hold listeners of events in the system
    this.markerClassName = 'hasEvent';
  }

  $.extend(Events.prototype, {
    // TODO rename to "widget" when switching to widget factory
    addSubscriber: function (name, subscriber) {
      //console.log('new subscriber [' + name + '] has been registered', subscriber);
      this.subscribers.push(subscriber);
    },

    /* This is used by boot upAttach the events to a jQuery selection.
     * @param  target	element - the target input field or division or span
     * @param  event event - The targets event we're interested in listening to
     * @param  settings  object - override settings to use on this event recording sequence (alternatives are data attrs)
     */
    attachEvents: function (el) {
      var $el = $(el);
      $el.addClass(this.markerClassName);
      if ($el.is('a')) {
        $el.click($.proxy(function (event) {
          var eventPayload = $el.data('cupEvent');
          if (eventPayload) {
            this.dispatchEvent(eventPayload);
          }
        }, this));
      }
    },
    /**
     * Dispatch an event to all current listeners of the events framework.
     *
     * @param ev - The event object, with an appropriate payload (see top of file).
     * @private
     */
    dispatchEvent: function (event, callback) {

      // this calls the callback should events not return in time
      var calledBack = false;
      var numberOfSubscribers = this.subscribers.length;
      var numberOfCallbacks = 0;

      function dispatch() {
        if (callback && !calledBack) {
          callback();
          calledBack = true;
        }
      }

      // start the race! -- this will get called if the subscribers don't callback in time!
      setTimeout(function () {
        dispatch(true);
      }, 3000);

      this.subscribers.forEach(function (subscriber) {
        subscriber.onEvent.call(this, event, function (err, res) {
          numberOfCallbacks++;
          if (numberOfCallbacks === numberOfSubscribers) {
            dispatch();
          }
        });
      });
    }
  });


  /* Invoke the events functionality.
   @param  options  string - a command, optionally followed by additional parameters or
   Object - settings for attaching new event functionality
   @return  jQuery object */
  $.fn.cupEvent = function () {
    $(this).find('[data-cup-event]').not('.hasEvent').each(function (i, el) {
      $.cupEvent.attachEvents(el);
    });
  };

  $.cupEvent = new Events(); // singleton instance
  return $.cupEvent;
})(jQuery);

