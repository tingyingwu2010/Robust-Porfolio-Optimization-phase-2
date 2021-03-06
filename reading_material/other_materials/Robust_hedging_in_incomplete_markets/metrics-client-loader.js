/* global google */
var AOP = AOP || {};

AOP.metrics = {

  $metricsTabContent: null,
  $metricsLoader: null,
  $metricsContentArea: null,
  $metricsTabButton: null,
  $metricsChartWrapper: null,
  $metricsChartEmpty: null,
  $metricsSummaryWrapper: null,

  chartsToRender: [],
  chartCount: 0,
  metricsLoaded: false,

  lookup: {
    service: AOP.baseUrl + '/services/aop-cambridge-core/metrics/',
    productId: null
  },


  init: function () {

    var _self = this;

    _self.limitSingleChartWidth = true;
    _self.$metricsTabContent = $('#metrics');
    _self.$metricsLoader = this.$metricsTabContent.find('.loader');
    _self.$metricsTabContentArea = this.$metricsTabContent.find('.reading-width');
    _self.$metricsTabButton = $('.tab-title a[href="#metrics"]');
    _self.$metricsChartWrapper = $('#metricsChartWrapper');
    _self.$metricsChartEmpty = $('#metricsChartEmpty');
    _self.$metricsSummaryWrapper = $('#metricsSummaryWrapper');
    _self.lookup.productId = this.$metricsTabContent.attr('data-product-id');

    _self.$metricsTabButton.on('click', function () {
      // google charts cannot render a graph correctly unless the element is visible in the DOM.
      // delayed loading of data now removes the need for setting a timeout
      if ( !_self.metricsLoaded) {
        _self.lookupData();
      }
    });
    if (_self.$metricsTabContent.is(':visible')) {
      _self.lookupData();
    }
  },


  /**
   * Retrieve data from metrics service endpoint
   */
  lookupData: function () {
    var _self = this;
    var width = this.getTargetChartWidth();
    var _query = '?_csrf=' + $('input[name="_csrf"]').val();
    var serviceURL = this.lookup.service + this.lookup.productId + _query;

    _self.showLoader(width);

    $.get(serviceURL).done(function (resp) {
      if (!resp.success) {
        _self.handleSadPath();
      } else {
        _self.chartsToRender = AOP.metrics.buildDataStructures(resp);
        _self.countChartsToRender();
        _self.renderCharts();
        _self.metricsLoaded = true;
      }
    });
  },


  /**
   * Build data into usable chart format
   */
  buildDataStructures: function(rawData) {

    //rawData = undefined;

    var rawFullText = {};
    var rawAbstract = {};

    if (rawData && rawData.data) {
      rawFullText = rawData.data.fullTextViews || {};
      rawAbstract = rawData.data.abstractViews || {};
    } else {
      this.handleSadPath();
      return [];
    }

    var rawFullTextData = rawFullText.data && rawFullText.data.allMonths || {};
    var rawAbstractData = rawAbstract.data && rawAbstract.data.allMonths || {};

    var fullTextHtmlTotal = rawFullText.viewTotal && rawFullText.viewTotal.html + 0 || 0;
    var fullTextPDFTotal = rawFullText.viewTotal && rawFullText.viewTotal.pdf + 0 || 0;


    var opts = {
      fullTextChartData: [],
      abstractChartData: [],
      fullTextTotal: parseInt(fullTextHtmlTotal) + parseInt(fullTextPDFTotal),
      abstractViewTotal: rawAbstract.viewTotal + 0
    };

    if (opts.fullTextTotal > 0 && rawFullTextData) {
      $('#fullTextViewsHTML').html(fullTextHtmlTotal);
      $('#fullTextViewsPDF').html(fullTextPDFTotal);

      for ( var month in rawFullTextData ) {
        opts.fullTextChartData.push([
          {v:month},
          rawFullTextData[month].totalAccumulated
        ]);
      }
    } else {
      $('#metricsFullTextWrapper').css('display', 'none');
    }

    if (opts.abstractViewTotal > 0 && rawAbstractData) {
      $('#abstractViewTotal').html(opts.abstractViewTotal);

      for ( var month in rawAbstractData ) {
        opts.abstractChartData.push([
          {v:month},
          rawAbstractData[month].viewsAccumulated
        ]);
      }
    } else {
      $('#metricsAbstractWrapper').css('display', 'none');
    }

    if (fullTextHtmlTotal === 0 && opts.abstractViewTotal === 0) {
      this.handleSadPath();
      return [];
    }

    this.upDateSummaryText(rawData);

    return this.createChartsToRender(opts);
  },

  handleSadPath: function () {
    this.$metricsChartWrapper.hide();
    this.$metricsChartEmpty.show();
    this.hideLoader();
  },

  upDateSummaryText: function (rawData){
    var dateRange = rawData.data && rawData.data.capturedDateRangeDescription || undefined;

    if (dateRange) {

      var $defaultSummaryElement = $('#defaultMetricSummary');
      var $multiLingSummaryElement = $('#multiLingualMetricSummary');

      if ($defaultSummaryElement.length) {
        var defaultSummaryText = $defaultSummaryElement.html();
        $defaultSummaryElement.html(defaultSummaryText.replace('#date#', dateRange));
      }

      if ($multiLingSummaryElement.length) {
        var multiLingDateRange = AOP.translateDateString(dateRange);
        var multiLingSummaryText = $multiLingSummaryElement.html();
        $multiLingSummaryElement.html(multiLingSummaryText.replace('&lt;date&gt;', multiLingDateRange));
      }
    } else {
      this.$metricsSummaryWrapper.hide();
    }
  },


  /**
   * Build chart properties
   */
  createChartsToRender: function(opts) {

    var hAxisLabel = this.$metricsChartWrapper.attr('data-charts-label-hAxis');
    var vAxisLabel = this.$metricsChartWrapper.attr('data-charts-label-vAxis');

    return [
      {
        containerId: 'article-months-chart-fulltext-views',
        chartType: 'LineChart',
        data: opts.fullTextChartData,
        total: opts.fullTextTotal,
        columns: [
          ['string', 'Month'],
          ['number', 'Full text views']
        ],
        chartOptions: {
          hAxis: {
            title: hAxisLabel
          },
          vAxis: {
            title: vAxisLabel,
            format: 0
          },
          legend: 'none',
          lineWidth: 3,
          pointSize: 5
        }
      },
      {
        containerId: 'article-months-chart-abstract-views',
        data: opts.abstractChartData,
        total: opts.abstractViewTotal,
        chartType: 'LineChart',
        columns: [
          ['string', 'Month'],
          ['number', 'Abstract views']
        ],
        chartOptions: {
          hAxis: {
            title: hAxisLabel
          },
          vAxis: {
            title: vAxisLabel,
            format: 0
          },
          legend: 'none',
          lineWidth: 3,
          pointSize: 5
        }
      }
    ]
  },


  /**
   * Count how many valid charts have been passed in
   */
  countChartsToRender: function () {
    var _self = this;
    _self.chartsToRender.forEach(function (chart) {
      if (chart.data && chart.data.length > 0) {
        _self.chartCount++;
      }
    });
  },


  /**
   * Set loader size and display
   */
  showLoader: function (screenWidth) {
    if (screenWidth < 200) {
      screenWidth = 200
    }
    this.$metricsLoader
      .css({width: screenWidth, height: Math.ceil(screenWidth / 2)})
      .find('>div').fadeIn();
  },


  /**
   * Hide the loader
   */
  hideLoader: function () {
    this.$metricsLoader.hide();
  },


  /**
   * Get the width the area to render the chart into
   * @returns {*|jQuery}
   */
  getTargetChartWidth: function () {
    // Get current chart container width
    var chartWidth = this.$metricsTabContentArea.width() || 300;

    // If we only have the single chart AND we are dealing with screen size larger than a standard table, reduce
    // the width the the graph to improve screen visuals.
    if (this.limitSingleChartWidth && this.chartCount === 1 && window.innerWidth > 1024) {
      chartWidth = chartWidth * 0.7;
    }
    return chartWidth;
  },


  /**
   * Render a single chart
   * @param opts
   */
  renderSingleChart: function (opts) {
    var _self = this;
    if (opts.data && opts.data.length > 0) {
      var chart = AOP.charts.initChart(opts);
      chart.draw();
      // If the user has re-sized the page while the tab has not been visible,
      // need re-render the charts so they fit correctly on the page.
      $(window).resize(function () {
        var chartWidth = _self.getTargetChartWidth();
        chart.draw({width: chartWidth});
      });
    }
  },


  /**
   * Render all the supplied charts
   */
  renderCharts: function () {
    if (this.chartsToRender.length) {
      var _self = this;
      var width = _self.getTargetChartWidth();
      AOP.charts.load({packages: ['corechart', 'bar']}, function () {
        _self.hideLoader();
        _self.chartsToRender.forEach(function (chart) {
          chart.width = width;
          _self.renderSingleChart(chart);
        });
      });
    }
  }
};

$(document).ready(AOP.metrics.init());
