/* globals AOP:true, Webtrends */

// WebTrends SmartSource Data Collector Tag v10.4.23
// Copyright (c) 2016 Webtrends Inc.  All rights reserved.
// Tag Builder Version: 4.1.3.5
// Created: 2016.01.13
window.webtrendsAsyncInit = function () {
  AOP = AOP || {};
  var dcs = new Webtrends.dcs().init({
    dcsid: AOP.webtrendsSourceId,
    domain: 'statse.webtrendslive.com',
    timezone: 0,
    i18n: true,
    download: true,
    downloadtypes: 'xls,doc,pdf,txt,csv,zip,docx,xlsx,rar,gzip',
    fpcdom: AOP.webTrendsFpcdom,
    plugins: {
      //hm:{src:'//s.webtrends.com/js/webtrends.hm.js'}
      omp: {
        src: AOP.baseUrl + '/cambridge-core/public/js/OracleMigrationPlugin.js',
        destinations: [
          {accountGuid: AOP.oracleInfinityAccountId, server: 'dc.oracleinfinity.io'}
        ],
        'waitForCallback': true,
        Qb: true
      }
    }
  }).track();
};
(function () {
  var s = document.createElement('script');
  s.async = true;
  s.src = AOP.baseUrl + '/cambridge-core/public/js/webtrends.min.js';

  var s2 = document.getElementsByTagName('script')[0];
  s2.parentNode.insertBefore(s, s2);
}());
