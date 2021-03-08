var FRJournalArticle=function(){function y(){$(window).on("hashchange",function(){o()});$(document).on("click",".ui-accordion",function(n){$(n.target).toggleClass("closed").next(".ui-accordion-content").slideToggle()});$(document).on("click","#supplementary_view",function(){s(!1);h()});$(document).on("click",function(n){if($(n.target).parents().andSelf().is(".popover, .popover-trigger")){var i=t.find(".popover");i.length>1&&i.each(function(){$(this).parent().has(n.target).length||$(this).prev(".popover-trigger").popover("hide")})}else $(t.find("a.popover-trigger")).popover("hide")});$(document).on("click",'a[href^="#"]',function(){var n=$(this).attr("href");n.length&&(n==="#impact"||n.indexOf("-tab")>=0||a(!0,n))});n.on("click","button[data-href]",function(){var n=$(this).data("href");n&&n.length&&(window.location.href=n)});n.find(".table-of-contents ul li a").on("click",function(){var n=$(this).attr("href");n&&n.length&&(window.location.href=n)});r.find("#open-crossmark").on("click",function(){});rt();nt()}function o(){var n=window.location.hash.replace("#","").toLowerCase();if(n.length&&n=="supplementary-material")h();else{f.modal("hide");return}}function p(){t.find(".side-article-impact .impact-data").tooltip({placement:"bottom",delay:{show:400}});$.getJSON("/articles/getarticleimpactviewscount?articleId="+window.FRArticle.ArticleId).done(function(n){t.find(".impact-data span.title-number").text(n)})}function w(){$.post("/articles/recommendedarticle/insertarticlerecommendations/"+window.FRArticle.ArticleId)}function b(){$.getJSON("/articles/getrelatedarticlebyarticleId?articleId="+window.FRArticle.ArticleId).done(function(n){k(n)})}function s(n,t){var i=t||1;setTimeout(function(){var t=$(".page-container a.open-supplemental-data");n?t.removeClass("disabled"):t.addClass("disabled")},i)}function h(){var n=FRTemplate.bind("template-supplementary-files-modal",{});f.html(n).modal("show");e=$("#localFiles");i=document.getElementById("loader");s(!0,800);i.style.display="block";$.ajax({type:"POST",url:FRArticle.FigShareApiUrl,data:JSON.stringify({doi:FRArticle.DOI}),dataType:"json",contentType:"application/json; charset=utf-8",timeout:FRArticle.FigShareTimeOut,success:function(n){n&&n.length!=0?ut():v()},error:function(){v()}})}function k(n){var t,i,r,f;n!=null&&(t=u.find(".widget-listing.original-article"),i=u.find(".widget-listing.commentary-article"),n.IsCommentaryVisible&&(r=c(n.CommentaryArticle),i.append(r).removeClass("hidden")),n.IsOriginalArticleVisible&&(f=c(n.OriginalArticle),t.append(f).removeClass("hidden")))}function c(n){var t="";return $.each(n,function(n,i){t=t+'<div class="teaser"><div class="heading-container"><h2 class="teaser-heading"><a href="'+i.ArticleUrl+'" >'+formatTextContent(i.ArticleTitle)+"<\/a><\/h2>";t=t+'<p class="teaser-authors">'+l(i.AuthorDetails)+"<\/p><\/div><\/div>"}),t}function d(){$.getJSON("/articles/recommendedarticle/getrecommendedarticledetails?articleId="+window.FRArticle.ArticleId).done(function(n){g(n)})}function g(n){if(n&&n.RecommendedArticleDetails&&n.RecommendedArticleDetails.length){var t=u.find(".widget-listing.people-also-looked-at");$.each(n.RecommendedArticleDetails,function(n,i){var r=$('<div class="teaser"><div class="heading-container"><h2 class="teaser-heading"><a><\/a><\/h2><p class="teaser-authors"><\/p><\/div><\/div>');r.find("a").attr("href",i.ArticleUrl).html(formatTextContent(i.ArticleTitle));r.find("p.teaser-authors").html(l(i.AuthorDetails));t.append(r)});t.removeClass("hidden")}}function l(n){if(!n||!n.length)return"";var i=0,t="",r=n.length;return $.each(n,function(n,u){(i++,t=u.WhosWhoUrl.length?t+'<a href="'+u.WhosWhoUrl+'">'+u.FullName+"<\/a>":t+u.FullName,r!==1)&&(i==r-1?t=t+" and ":i!=r&&(t=t+", "))}),t}function nt(){var n=t.find("div.altmetric-icon");n.on("click",".altmetric-embed>a",function(){})}function tt(){r.find('a[id^="h"]').addClass("reset-hash-position")}function a(n,t){var i,u,f;t&&t.length&&(i=t.indexOf("/"),i>0&&(t=t.substring(0,i)),u=function(){var n=r.find(t),i;n&&n.length&&(n.hasClass("reset-hash-position")||(i=parseInt(n.offset().top)-parseInt($("body").css("padding-top")),$("html, body").animate({scrollTop:i},100)))},f=n?50:500,setTimeout(u,f))}function it(){n.find("#top-impact-factor").circleProgress({value:1,size:50,thickness:5,fill:{color:"rgba(0, 0, 0, 0.3)"}}).on("circle-animation-progress",function(){});n.find("#top-impact-factor .counter").counterUp()}function rt(){var r;n.find("#anchors").find("a[href*=#]:not([href=#], .nav-tabs a)").click(function(){if(location.pathname.replace(/^\//,"")!=this.pathname.replace(/^\//,"")||location.hostname!=this.hostname)return!1;var n=$(this.hash);return(n=n.length?n:$("[name="+this.hash.slice(1)+"]"),!n.length)?!1:($("html,body").animate({scrollTop:n.offset().top},1e3),!1)});n.find("#anchors").stick_in_parent();var t=200,i=!1,u=new Date(1,1,2e3,12,00,00);$(window).resize(function(){u=new Date;i===!1&&(i=!0,setTimeout(r,t))});r=function(){new Date-u<t?setTimeout(r,t):(i=!1,$(document.body).trigger("sticky_kit:recalc"))}}function ut(){window.figshare.load("frontiers",function(n){var r=document.getElementById("figshare-widget-container"),t=new n({resourceDOI:window.FRArticle.DOI,showStats:!1});t.initialize();t.mount(r);window.widget=t;i.style.display="none"})}function v(){$.getJSON("/articles/GetSupplementaryFilesByArticleId?articleId="+window.FRArticle.ArticleId).done(function(n){var t=n.SupplimentalFileDetails,r=FRTemplate.bind("template-local-files-modal",{supplimentalFileDetails:t});e.html(r);i.style.display="none";t.FileDetails==null?($(".supplementary-content").hide(),$(".supplementary-empty-message").show()):($(".supplementary-content").show(),$(".supplementary-empty-message").hide())})}var n=$(".page-container"),t=n.find(".right-container"),r=n.find(".article-container"),ft=n.find(".abstract-container"),u=n.find(".right-container-articles"),f=$(".supplementary-modal-container"),e=null,i=null;$(function(){p();b();d();w();y();tt();a(!1,window.location.hash);FRNetworkUserFollow.fillJournalFollow();it();o()})}(),FRAddThis=function(){function h(){window.addthis&&window.addthis.addEventListener&&window.addthis.addEventListener("addthis.menu.share",c)}function c(n){return n.type=="addthis.menu.share"&&p(n.data.service),!0}function l(){var n=t.find(".addthis_button_facebook");i.on("click",function(){if(n.length)n.click();else return});u.on("click",function(){var n=t.find(".addthis_button_google");if(n.length)n.click();else return});r.on("click",function(){var n=t.find(".addthis_button_twitter");if(n.length)n.click();else return});f.on("click",function(){var n=t.find(".addthis_button_linkedin");if(n.length)n.click();else return});e.on("click",function(){var n=t.find(".addthis_button_more");if(n.length)n.click();else return})}function a(){window.FRSocial.pageType==1?v():window.FRSocial.pageType==2&&y()}function v(){$.getJSON("/articles/social/getsocialcountbyarticleId?articleid="+window.FRArticle.ArticleId).done(function(n){s(n)})}function y(){$.getJSON("/research-topics/"+window.FRResearchTopic.Id+"/socialcounts").done(function(n){s(n)})}function s(n){typeof n!="undefined"&&n!=null&&(typeof n.Facebook=="undefined"?(i.html(0),i.data("count",0)):(i.html(n.Facebook),i.data("count",n.Facebook)),typeof n.Twitter=="undefined"?(r.html(0),r.data("count",0)):(r.html(n.Twitter),r.data("count",n.Twitter)),typeof n.GooglePlus=="undefined"?(u.html(0),u.data("count",0)):(u.html(n.GooglePlus),u.data("count",n.GooglePlus)),typeof n.LinkedIn=="undefined"?(f.html(0),f.data("count",0)):(f.html(n.LinkedIn),f.data("count",n.LinkedIn)),typeof n.Total=="undefined"||n.Total==0?(e.html("New"),e.data("count",0)):(e.html(n.Total),e.data("count",n.Total)))}function p(t){var i;t!="more"&&(t=="facebook"?(i=n.find(".facebook_count"),o(i)):t=="twitter"?(i=n.find(".twitter_count"),o(i)):t=="google"||t=="google_plusone_share"?(t="google",i=n.find(".googleplus_count"),o(i)):t=="linkedin"?(i=n.find(".linkedin_count"),o(i)):k(),window.FRSocial.pageType==1?w(t):window.FRSocial.pageType==2&&b(t))}function w(n){$.post("/Social/TrackCurrentSocialCounts?id="+window.FRArticle.ArticleId+"&provider="+n)}function b(n){$.post("/ResearchTopic/TrackCurrentSocialCounts?id="+window.FRResearchTopic.Id+"&provider="+n)}function o(t){var u=$(t).data("count"),i=parseInt(u)+parseInt(1),f=n.find(".total_count").data("count"),r=parseInt(f)+parseInt(1);$(t).html(i);n.find(".total_count").html(r);$(t).data("count",i);n.find(".total_count").data("count",r)}function k(){var i=n.find(".total_count").data("count"),t=parseInt(i)+parseInt(1);n.find(".total_count").html(t);n.find(".total_count").data("count",t)}var n=$(".share-media"),t=n.find(".social_block"),i=t.find(".facebook_count"),r=t.find(".twitter_count"),u=t.find(".googleplus_count"),f=t.find(".linkedin_count"),e=t.find(".total_count");$(function(){h();a();l()})}(),onloadCallback=function(){grecaptcha.render("recaptcha",{sitekey:window.FRArticleRecaptchaSettings.RecaptchaSiteKey,callback:function(){FRArticleNotification.SetNotifyVisibility()},"expired-callback":function(){FRArticleNotification.SetNotifyVisibility()}})},FRArticleNotification=function(){function e(){(function(n,t,i,r,u,f,e){f=t.createElement(i);e=t.getElementsByTagName(i)[0];f.async=1;f.defer=1;f.src=r;e.parentNode.insertBefore(f,e)})(window,document,"script","https://www.google.com/recaptcha/api.js?onload=onloadCallback&render=explicit")}function o(i){var r=i.currentTarget;r.innerHTML="";$.get("/api/articles/getloggedinuser").done(function(i){if(i>0)$.get("/api/articles/"+i+"/getuseremaillist").done(function(r){if(r.length>1){var u=FRTemplate.bind("template-notifyme-mailselection-modal",{userEmailList:r});n.html(u).modal("show");$('input[name="notify_Emailcollection"]').first().prop("checked",!0)}else $.post("/api/articles/notifyme",{ArticleId:window.FRArticle.ArticleId,EmailId:r[0],UserId:i}).done(function(n){n?$(".provisional-text").text("We'll notify you at publication."):t()})});else{var r=FRTemplate.bind("template-notify-me-modal",{login:FRConfiguration,currentPage:location.pathname});n.html(r).modal("show");e()}})}function t(t){t=t||f;var i=FRTemplate.bind("template_notifyme_error_modal",{message:t});n.modal("hide");setTimeout(function(){n.html(i).modal("show")},500)}function i(n,i,u){var f=n.currentTarget;f.innerHTML='<i class="fa fa-spin fa-circle-o-notch"><\/i>';u?$.post("/api/articles/notifyme/",{ArticleId:window.FRArticle.ArticleId,EmailId:i,UserId:FRArticle.LoginUserId}).done(function(n){r(n,f)}).fail(function(){t()}):$.post("/api/articles/notification/security/verification/",{ArticleId:window.FRArticle.ArticleId,EmailId:i,UserId:FRArticle.LoginUserId,RecaptchaToken:grecaptcha.getResponse()}).done(function(n){r(n,f)}).fail(function(){t("Your Captcha response was incorrect. Please try again.")})}function r(i,r){i?(r.innerHTML='<i class="fa fa-check" aria-hidden="true"><\/i>',n.modal("hide"),FRArticle.LoginUserId>0&&$(".provisional-text").text("We'll notify you at publication.")):t()}function s(){var n=!1;return typeof grecaptcha!="undefined"&&grecaptcha.getResponse()&&(n=!0),n}function u(){var n=s(),t=n?"visible":"hidden";$("#article_notify_non_registered_user").css("visibility",t)}var n=$(".notifyme-modal-container"),h=$("#article_notify_loggedin_user"),f="An error occured. Please try again later.";$(function(){u()});$(document).on("click","#article_notifyme",function(n){o(n)});$(document).on("click","#article_notify_non_registered_user",function(n){if($("#modalnotifyme").valid()){var t=$("#txt_notification_email_id").val();i(n,t,!1)}});$(document).on("click","#article_notify_loggedin_user",function(n){email=$("input[name='notify_Emailcollection']:checked").val();i(n,email,!0)});n.on("hidden.bs.modal",function(){$("#article_notifyme").html('<i class="fa fa-envelope"><\/i> <b>Notify me<\/b>')});return{SetNotifyVisibility:u}}(),JournalBranding=$(document).ready(function(){function n(){var n=$("#journal-header .associate-inner img.associate-inner__image");n&&t(n[0])}function t(n){let t="associate-inner__image";if(n){var i=n?n.naturalHeight:0,r=n?n.naturalWidth:0;i>=r&&(t=t+" portrait-orientation");n.setAttribute("class",t)}}return $(function(){n()}),{}});(function(n,t){"object"==typeof exports&&"object"==typeof module?module.exports=t():"function"==typeof define&&define.amd?define([],t):"object"==typeof exports?exports.footer=t():n.footer=t()})(window,function(){var n=Math.min;return function(n){function t(r){if(i[r])return i[r].exports;var u=i[r]={i:r,l:!1,exports:{}};return n[r].call(u.exports,u,u.exports,t),u.l=!0,u.exports}var i={};return t.m=n,t.c=i,t.d=function(n,i,r){t.o(n,i)||Object.defineProperty(n,i,{enumerable:!0,get:r})},t.r=function(n){"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(n,Symbol.toStringTag,{value:"Module"});Object.defineProperty(n,"__esModule",{value:!0})},t.t=function(n,i){var r,u;if((1&i&&(n=t(n)),8&i)||4&i&&"object"==typeof n&&n&&n.__esModule)return n;if(r=Object.create(null),t.r(r),Object.defineProperty(r,"default",{enumerable:!0,value:n}),2&i&&"string"!=typeof n)for(u in n)t.d(r,u,function(t){return n[t]}.bind(null,u));return r},t.n=function(n){var i=n&&n.__esModule?function(){return n["default"]}:function(){return n};return t.d(i,"a",i),i},t.o=function(n,t){return Object.prototype.hasOwnProperty.call(n,t)},t.p="",t(t.s=31)}([function(n,t,i){(function(t){var i=function(n){return n&&n.Math==Math&&n};n.exports=i("object"==typeof globalThis&&globalThis)||i("object"==typeof window&&window)||i("object"==typeof self&&self)||i("object"==typeof t&&t)||Function("return this")()}).call(this,i(34))},function(n){n.exports=function(n){try{return!!n()}catch(n){return!0}}},function(n){n.exports=function(n){return"object"==typeof n?null!==n:"function"==typeof n}},function(n){var t={}.hasOwnProperty;n.exports=function(n,i){return t.call(n,i)}},function(n,t,i){var e=i(0),o=i(24),f=i(3),s=i(25),h=i(29),c=i(55),u=o("wks"),r=e.Symbol,l=c?r:r&&r.withoutSetter||s;n.exports=function(n){return f(u,n)||(u[n]=h&&f(r,n)?r[n]:l("Symbol."+n)),u[n]}},function(n,t,i){var r=i(7),u=i(14),f=i(11);n.exports=r?function(n,t,i){return u.f(n,t,f(1,i))}:function(n,t,i){return n[t]=i,n}},function(n,t,i){var r=i(2);n.exports=function(n){if(!r(n))throw TypeError(n+" is not an object");return n}},function(n,t,i){var r=i(1);n.exports=!r(function(){return 7!=Object.defineProperty({},1,{get:function(){return 7}})[1]})},function(n){var t={}.toString;n.exports=function(n){return t.call(n).slice(8,-1)}},function(n){n.exports=function(n){if(n==void 0)throw TypeError("Can't call method on "+n);return n}},function(n,t,i){"use strict";var h=i(61),e=i(62),r=RegExp.prototype.exec,c=String.prototype.replace,o=r,u=function(){var n=/a/,t=/b*/g;return r.call(n,"a"),r.call(t,"a"),0!==n.lastIndex||0!==t.lastIndex}(),s=e.UNSUPPORTED_Y||e.BROKEN_CARET,f=/()??/.exec("")[1]!==void 0;(u||f||s)&&(o=function(n){var w,l,t,o,i=this,y=s&&i.sticky,e=h.call(i),a=i.source,p=0,v=n;return y&&(e=e.replace("y",""),-1===e.indexOf("g")&&(e+="g"),v=(n+"").slice(i.lastIndex),0<i.lastIndex&&(!i.multiline||i.multiline&&"\n"!==n[i.lastIndex-1])&&(a="(?: "+a+")",v=" "+v,p++),l=new RegExp("^(?:"+a+")",e)),f&&(l=new RegExp("^"+a+"$(?!\\s)",e)),u&&(w=i.lastIndex),t=r.call(y?l:i,v),y?t?(t.input=t.input.slice(p),t[0]=t[0].slice(p),t.index=i.lastIndex,i.lastIndex+=t[0].length):i.lastIndex=0:u&&t&&(i.lastIndex=i.global?t.index+t[0].length:w),f&&t&&1<t.length&&c.call(t[0],l,function(){for(o=1;o<arguments.length-2;o++)void 0===arguments[o]&&(t[o]=void 0)}),t});n.exports=o},function(n){n.exports=function(n,t){return{enumerable:!(1&n),configurable:!(2&n),writable:!(4&n),value:t}}},function(n,t,i){var r=i(36),u=i(9);n.exports=function(n){return r(u(n))}},function(n,t,i){var r=i(2);n.exports=function(n,t){if(!r(n))return n;var i,u;if(t&&"function"==typeof(i=n.toString)&&!r(u=i.call(n))||"function"==typeof(i=n.valueOf)&&!r(u=i.call(n))||!t&&"function"==typeof(i=n.toString)&&!r(u=i.call(n)))return u;throw TypeError("Can't convert object to primitive value");}},function(n,t,i){var f=i(7),e=i(20),r=i(6),o=i(13),u=Object.defineProperty;t.f=f?u:function(n,t,i){if(r(n),t=o(t,!0),r(i),e)try{return u(n,t,i)}catch(n){}if("get"in i||"set"in i)throw TypeError("Accessors not supported");return"value"in i&&(n[t]=i.value),n}},function(n,t,i){var r=i(0),u=i(5);n.exports=function(n,t){try{u(r,n,t)}catch(i){r[n]=t}return t}},function(t,i,r){var u=r(17);t.exports=function(t){return 0<t?n(u(t),9007199254740991):0}},function(n){var t=Math.ceil,i=Math.floor;n.exports=function(n){return isNaN(n=+n)?0:(0<n?i:t)(n)}},function(n,t,i){var r=i(0),u=i(19).f,f=i(5),e=i(21),o=i(15),s=i(42),h=i(51);n.exports=function(n,t){var p,l,i,c,a,y,v=n.target,w=n.global,b=n.stat;if(l=w?r:b?r[v]||o(v,{}):(r[v]||{}).prototype,l)for(i in t){if(a=t[i],n.noTargetGet?(y=u(l,i),c=y&&y.value):c=l[i],p=h(w?i:v+(b?".":"#")+i,n.forced),!p&&void 0!==c){if(typeof a==typeof c)continue;s(a,c)}(n.sham||c&&c.sham)&&f(a,"sham",!0);e(l,i,a,n)}}},function(n,t,i){var u=i(7),f=i(35),e=i(11),o=i(12),s=i(13),h=i(3),c=i(20),r=Object.getOwnPropertyDescriptor;t.f=u?r:function(n,t){if(n=o(n),t=s(t,!0),c)try{return r(n,t)}catch(n){}if(h(n,t))return e(!f.f.call(n,t),n[t])}},function(n,t,i){var r=i(7),u=i(1),f=i(37);n.exports=!r&&!u(function(){return 7!=Object.defineProperty(f("div"),"a",{get:function(){return 7}}).a})},function(n,t,i){var f=i(0),r=i(5),e=i(3),o=i(15),s=i(22),u=i(38),h=u.get,c=u.enforce,l=(String+"").split("String");(n.exports=function(n,t,i,u){var h=!!u&&!!u.unsafe,s=!!u&&!!u.enumerable,a=!!u&&!!u.noTargetGet;return("function"==typeof i&&("string"==typeof t&&!e(i,"name")&&r(i,"name",t),c(i).source=l.join("string"==typeof t?t:"")),n===f)?void(s?n[t]=i:o(t,i)):void(h?!a&&n[t]&&(s=!0):delete n[t],s?n[t]=i:r(n,t,i))})(Function.prototype,"toString",function(){return"function"==typeof this&&h(this).source||s(this)})},function(n,t,i){var r=i(23),u=Function.toString;"function"!=typeof r.inspectSource&&(r.inspectSource=function(n){return u.call(n)});n.exports=r.inspectSource},function(n,t,i){var u=i(0),f=i(15),r="__core-js_shared__",e=u[r]||f(r,{});n.exports=e},function(n,t,i){var u=i(41),r=i(23);(n.exports=function(n,t){return r[n]||(r[n]=t===void 0?{}:t)})("versions",[]).push({version:"3.6.2",mode:u?"pure":"global",copyright:"© 2020 Denis Pushkarev (zloirock.ru)"})},function(n){var t=0,i=Math.random();n.exports=function(n){return"Symbol("+((n===void 0?"":n)+"")+")_"+(++t+i).toString(36)}},function(n){n.exports={}},function(n,t,i){var r=i(44),u=i(0),f=function(n){if("function"==typeof n)return n};n.exports=function(n,t){return 2>arguments.length?f(r[n])||f(u[n]):r[n]&&r[n][t]||u[n]&&u[n][t]}},function(n,t,i){var r=i(8);n.exports=Array.isArray||function(n){return"Array"==r(n)}},function(n,t,i){var r=i(1);n.exports=!!Object.getOwnPropertySymbols&&!r(function(){return!(Symbol()+"")})},function(n,t,i){var r,u,h=i(0),f=i(57),e=h.process,o=e&&e.versions,s=o&&o.v8;s?(r=s.split("."),u=r[0]+r[1]):f&&(r=f.match(/Edge\/(\d+)/),(!r||74<=r[1])&&(r=f.match(/Chrome\/(\d+)/),r&&(u=r[1])));n.exports=u&&+u},function(n,t,i){n.exports=i(32)},function(n,t,i){"use strict";i.r(t);var r=i(33),f=i.n(r),u=i(58),e=i.n(u);t["default"]={Render:function(n,t){function r(n,t){var r=document.getElementsByTagName("body")[0],i=document.createElement("script");i.src=n;t&&(i.onload=function(){return t(i)});r.appendChild(i)}function f(n){var i=document.getElementsByTagName("head")[0],t=document.createElement("link");t.rel="stylesheet";t.href=n;i.appendChild(t)}var u="",i;if(!t){if(i=window.location.hostname.split(".").reverse(),2>i.length)throw"Invalid hostname ".concat(window.location.hostname);u="".concat(i[1],".").concat(i[0]);t={footerApiUrl:"https://navigation-plugins-api.".concat(u,"/plugins-api/get-footer-resource")}}r(t.footerApiUrl,function(){f(footerNavigationPluginResource.info.StyleUrl);r(footerNavigationPluginResource.info.ScriptUrl,function(){footerNavigationPlugin.Render(n,footerNavigationPluginResource.menuConfiguration)})})}}},function(n,t,i){"use strict";var o=i(18),s=i(1),h=i(28),c=i(2),l=i(52),a=i(16),r=i(53),v=i(54),y=i(56),p=i(4),w=i(30),u=p("isConcatSpreadable"),f=9007199254740991,e="Maximum allowed index exceeded",b=51<=w||!s(function(){var n=[];return n[u]=!1,n.concat()[0]!==n}),k=y("concat"),d=function(n){if(!c(n))return!1;var t=n[u];return t===void 0?h(n):!!t};o({target:"Array",proto:!0,forced:!b||!k},{concat:function(){for(var u,s,n,c=l(this),o=v(c,0),t=0,i=-1,h=arguments.length;i<h;i++)if(n=-1===i?c:arguments[i],d(n)){if(s=a(n.length),t+s>f)throw TypeError(e);for(u=0;u<s;u++,t++)u in n&&r(o,t,n[u])}else{if(t>=f)throw TypeError(e);r(o,t++,n)}return o.length=t,o}})},function(n){var t=function(){return this}();try{t=t||new Function("return this")()}catch(n){"object"==typeof window&&(t=window)}n.exports=t},function(n,t){"use strict";var i={}.propertyIsEnumerable,r=Object.getOwnPropertyDescriptor,u=r&&!i.call({1:2},1);t.f=u?function(n){var t=r(this,n);return!!t&&t.enumerable}:i},function(n,t,i){var r=i(1),u=i(8),f="".split;n.exports=r(function(){return!Object("z").propertyIsEnumerable(0)})?function(n){return"String"==u(n)?f.call(n,""):Object(n)}:Object},function(n,t,i){var f=i(0),u=i(2),r=f.document,e=u(r)&&u(r.createElement);n.exports=function(n){return e?r.createElement(n):{}}},function(n,t,i){var e,f,o,h=i(39),c=i(0),l=i(2),a=i(5),s=i(3),v=i(40),y=i(26),p=c.WeakMap,w=function(n){return o(n)?f(n):e(n,{})},u;if(h){var r=new p,b=r.get,k=r.has,d=r.set;e=function(n,t){return d.call(r,n,t),t};f=function(n){return b.call(r,n)||{}};o=function(n){return k.call(r,n)}}else u=v("state"),y[u]=!0,e=function(n,t){return a(n,u,t),t},f=function(n){return s(n,u)?n[u]:{}},o=function(n){return s(n,u)};n.exports={set:e,get:f,has:o,enforce:w,getterFor:function(n){return function(t){var i;if(!l(t)||(i=f(t)).type!==n)throw TypeError("Incompatible receiver, "+n+" required");return i}}}},function(n,t,i){var u=i(0),f=i(22),r=u.WeakMap;n.exports="function"==typeof r&&/native code/.test(f(r))},function(n,t,i){var u=i(24),f=i(25),r=u("keys");n.exports=function(n){return r[n]||(r[n]=f(n))}},function(n){n.exports=!1},function(n,t,i){var r=i(3),u=i(43),f=i(19),e=i(14);n.exports=function(n,t){for(var i,s=u(t),h=e.f,c=f.f,o=0;o<s.length;o++)i=s[o],r(n,i)||h(n,i,c(t,i))}},function(n,t,i){var r=i(27),u=i(45),f=i(50),e=i(6);n.exports=r("Reflect","ownKeys")||function(n){var t=u.f(e(n)),i=f.f;return i?t.concat(i(n)):t}},function(n,t,i){var r=i(0);n.exports=r},function(n,t,i){var r=i(46),u=i(49),f=u.concat("length","prototype");t.f=Object.getOwnPropertyNames||function(n){return r(n,f)}},function(n,t,i){var r=i(3),u=i(12),f=i(47).indexOf,e=i(26);n.exports=function(n,t){var i,s=u(n),h=0,o=[];for(i in s)!r(e,i)&&r(s,i)&&o.push(i);for(;t.length>h;)r(s,i=t[h++])&&(~f(o,i)||o.push(i));return o}},function(n,t,i){var u=i(12),f=i(16),e=i(48),r=function(n){return function(t,i,r){var h,s=u(t),c=f(s.length),o=e(r,c);if(n&&i!=i){for(;c>o;)if(h=s[o++],h!=h)return!0}else for(;c>o;o++)if((n||o in s)&&s[o]===i)return n||o||0;return!n&&-1}};n.exports={includes:r(!0),indexOf:r(!1)}},function(t,i,r){var u=r(17),f=Math.max;t.exports=function(t,i){var r=u(t);return 0>r?f(r+i,0):n(r,i)}},function(n){n.exports=["constructor","hasOwnProperty","isPrototypeOf","propertyIsEnumerable","toLocaleString","toString","valueOf"]},function(n,t){t.f=Object.getOwnPropertySymbols},function(n,t,i){var u=i(1),f=/#|\.prototype\./,r=function(n,t){var i=o[e(n)];return!(i!=h)||i!=s&&("function"==typeof t?u(t):!!t)},e=r.normalize=function(n){return(n+"").replace(f,".").toLowerCase()},o=r.data={},s=r.NATIVE="N",h=r.POLYFILL="P";n.exports=r},function(n,t,i){var r=i(9);n.exports=function(n){return Object(r(n))}},function(n,t,i){"use strict";var r=i(13),u=i(14),f=i(11);n.exports=function(n,t,i){var e=r(t);e in n?u.f(n,e,f(0,i)):n[e]=i}},function(n,t,i){var u=i(2),r=i(28),f=i(4),e=f("species");n.exports=function(n,t){var i;return r(n)&&(i=n.constructor,"function"==typeof i&&(i===Array||r(i.prototype))?i=void 0:u(i)&&(i=i[e],null===i&&(i=void 0))),new(void 0===i?Array:i)(0===t?0:t)}},function(n,t,i){var r=i(29);n.exports=r&&!Symbol.sham&&"symbol"==typeof Symbol.iterator},function(n,t,i){var r=i(1),u=i(4),f=i(30),e=u("species");n.exports=function(n){return 51<=f||!r(function(){var t=[],i=t.constructor={};return i[e]=function(){return{foo:1}},1!==t[n](Boolean).foo})}},function(n,t,i){var r=i(27);n.exports=r("navigator","userAgent")||""},function(t,i,r){"use strict";var s=r(59),h=r(63),c=r(6),e=r(9),l=r(64),a=r(66),v=r(16),o=r(68),y=r(10),p=r(1),w=[].push,f=4294967295,u=!p(function(){return!RegExp(f,"y")});s("split",2,function(t,i,r){var s;return s="c"=="abbc".split(/(b)*/)[1]||4!="test".split(/(?:)/,-1).length||2!="ab".split(/(?:ab)*/).length||4!=".".split(/(.?)(.?)/).length||1<".".split(/()()/).length||"".split(/.?/).length?function(n,t){var u=e(this)+"",s=void 0===t?f:t>>>0;if(0==s)return[];if(void 0===n)return[u];if(!h(n))return i.call(u,n,s);for(var o,a,v,r=[],p=(n.ignoreCase?"i":"")+(n.multiline?"m":"")+(n.unicode?"u":"")+(n.sticky?"y":""),c=0,l=new RegExp(n.source,p+"g");(o=y.call(l,u))&&(a=l.lastIndex,!(a>c&&(r.push(u.slice(c,o.index)),1<o.length&&o.index<u.length&&w.apply(r,o.slice(1)),v=o[0].length,c=a,r.length>=s)));)l.lastIndex===o.index&&l.lastIndex++;return c===u.length?(v||!l.test(""))&&r.push(""):r.push(u.slice(c)),r.length>s?r.slice(0,s):r}:function(n,t){return void 0===n&&0===t?[]:i.call(this,n,t)},[function(n,i){var r=e(this),u=void 0==n?void 0:n[t];return void 0===u?s.call(r+"",n,i):u.call(n,r,i)},function(t,e){var tt=r(s,t,this,e,s!==i),it,d,g;if(tt.done)return tt.value;var w=c(t),h=this+"",rt=l(w,RegExp),ut=w.unicode,ft=(w.ignoreCase?"i":"")+(w.multiline?"m":"")+(w.unicode?"u":"")+(u?"y":"g"),b=new rt(u?w:"^(?:"+w.source+")",ft),nt=void 0===e?f:e>>>0;if(0==nt)return[];if(0===h.length)return null===o(b,h)?[h]:[];for(var k=0,y=0,p=[];y<h.length;)if(b.lastIndex=u?y:0,d=o(b,u?h:h.slice(y)),null===d||(it=n(v(b.lastIndex+(u?0:y)),h.length))===k)y=a(h,y,ut);else{if(p.push(h.slice(k,y)),p.length===nt)return p;for(g=1;g<=d.length-1;g++)if(p.push(d[g]),p.length===nt)return p;y=k=it}return p.push(h.slice(k)),p}]},!u)},function(n,t,i){"use strict";i(60);var u=i(21),r=i(1),f=i(4),o=i(10),s=i(5),h=f("species"),c=!r(function(){var n=/./;return n.exec=function(){var n=[];return n.groups={a:"7"},n},"7"!=="".replace(n,"$<a>")}),e=function(){return"$0"==="a".replace(/./,"$0")}(),l=!r(function(){var t=/(?:)/,i=t.exec,n;return t.exec=function(){return i.apply(this,arguments)},n="ab".split(t),2!==n.length||"a"!==n[0]||"b"!==n[1]});n.exports=function(n,t,i,a){var v=f(n),y=!r(function(){var t={};return t[v]=function(){return 7},7!=""[n](t)}),b=y&&!r(function(){var i=!1,t=/a/;return"split"===n&&(t={},t.constructor={},t.constructor[h]=function(){return t},t.flags="",t[v]=/./[v]),t.exec=function(){return i=!0,null},t[v](""),!i});if(!y||!b||"replace"===n&&!(c&&e)||"split"===n&&!l){var k=/./[v],p=i(v,""[n],function(n,t,i,r,u){return t.exec===o?y&&!u?{done:!0,value:k.call(t,i,r)}:{done:!0,value:n.call(i,t,r)}:{done:!1}},{REPLACE_KEEPS_$0:e}),d=p[0],w=p[1];u(String.prototype,n,d);u(RegExp.prototype,v,2==t?function(n,t){return w.call(n,this,t)}:function(n){return w.call(n,this)})}a&&s(RegExp.prototype[v],"sham",!0)}},function(n,t,i){"use strict";var u=i(18),r=i(10);u({target:"RegExp",proto:!0,forced:/./.exec!==r},{exec:r})},function(n,t,i){"use strict";var r=i(6);n.exports=function(){var t=r(this),n="";return t.global&&(n+="g"),t.ignoreCase&&(n+="i"),t.multiline&&(n+="m"),t.dotAll&&(n+="s"),t.unicode&&(n+="u"),t.sticky&&(n+="y"),n}},function(n,t,i){"use strict";function r(n,t){return RegExp(n,t)}var u=i(1);t.UNSUPPORTED_Y=u(function(){var n=r("a","y");return n.lastIndex=2,null!=n.exec("abcd")});t.BROKEN_CARET=u(function(){var n=r("^r","gy");return n.lastIndex=2,null!=n.exec("str")})},function(n,t,i){var r=i(2),u=i(8),f=i(4),e=f("match");n.exports=function(n){var t;return r(n)&&((t=n[e])===void 0?"RegExp"==u(n):!!t)}},function(n,t,i){var r=i(6),u=i(65),f=i(4),e=f("species");n.exports=function(n,t){var i,f=r(n).constructor;return f===void 0||(i=r(f)[e])==void 0?t:u(i)}},function(n){n.exports=function(n){if("function"!=typeof n)throw TypeError(n+" is not a function");return n}},function(n,t,i){"use strict";var r=i(67).charAt;n.exports=function(n,t,i){return t+(i?r(n,t).length:1)}},function(n,t,i){var u=i(17),f=i(9),r=function(n){return function(t,i){var e,s,o=f(t)+"",r=u(i),h=o.length;return 0>r||r>=h?n?"":void 0:(e=o.charCodeAt(r),55296>e||56319<e||r+1===h||56320>(s=o.charCodeAt(r+1))||57343<s?n?o.charAt(r):e:n?o.slice(r,r+2):(e-55296<<10)+(s-56320)+65536)}};n.exports={codeAt:r(!1),charAt:r(!0)}},function(n,t,i){var r=i(8),u=i(10);n.exports=function(n,t){var f=n.exec,i;if("function"==typeof f){if(i=f.call(n,t),"object"!=typeof i)throw TypeError("RegExp exec method returned something other than an Object or null");return i}if("RegExp"!==r(n))throw TypeError("RegExp#exec called on incompatible receiver");return u.call(n,t)}}])["default"]}),function(n,t){"object"==typeof exports&&"object"==typeof module?module.exports=t():"function"==typeof define&&define.amd?define([],t):"object"==typeof exports?exports.header=t():n.header=t()}(window,function(){return function(n){function t(r){if(i[r])return i[r].exports;var u=i[r]={i:r,l:!1,exports:{}};return n[r].call(u.exports,u,u.exports,t),u.l=!0,u.exports}var i={};return t.m=n,t.c=i,t.d=function(n,i,r){t.o(n,i)||Object.defineProperty(n,i,{enumerable:!0,get:r})},t.r=function(n){"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(n,Symbol.toStringTag,{value:"Module"});Object.defineProperty(n,"__esModule",{value:!0})},t.t=function(n,i){var r,u;if((1&i&&(n=t(n)),8&i)||4&i&&"object"==typeof n&&n&&n.__esModule)return n;if(r=Object.create(null),t.r(r),Object.defineProperty(r,"default",{enumerable:!0,value:n}),2&i&&"string"!=typeof n)for(u in n)t.d(r,u,function(t){return n[t]}.bind(null,u));return r},t.n=function(n){var i=n&&n.__esModule?function(){return n["default"]}:function(){return n};return t.d(i,"a",i),i},t.o=function(n,t){return Object.prototype.hasOwnProperty.call(n,t)},t.p="",t(t.s=21)}([function(n,t,i){(function(t){var i=function(n){return n&&n.Math==Math&&n};n.exports=i("object"==typeof globalThis&&globalThis)||i("object"==typeof window&&window)||i("object"==typeof self&&self)||i("object"==typeof t&&t)||Function("return this")()}).call(this,i(27))},function(n){n.exports=function(n){try{return!!n()}catch(n){return!0}}},function(n){n.exports=function(n){return"object"==typeof n?null!==n:"function"==typeof n}},function(n,t,i){var r=i(1);n.exports=!r(function(){return 7!=Object.defineProperty({},1,{get:function(){return 7}})[1]})},function(n){n.exports=function(n,t){return{enumerable:!(1&n),configurable:!(2&n),writable:!(4&n),value:t}}},function(n,t,i){var r=i(2);n.exports=function(n,t){if(!r(n))return n;var i,u;if(t&&"function"==typeof(i=n.toString)&&!r(u=i.call(n))||"function"==typeof(i=n.valueOf)&&!r(u=i.call(n))||!t&&"function"==typeof(i=n.toString)&&!r(u=i.call(n)))return u;throw TypeError("Can't convert object to primitive value");}},function(n){var t={}.hasOwnProperty;n.exports=function(n,i){return t.call(n,i)}},function(n){n.exports={}},function(n,t,i){var r=i(11);n.exports=Array.isArray||function(n){return"Array"==r(n)}},function(n,t,i){var e=i(0),o=i(42),f=i(6),s=i(46),h=i(16),c=i(47),u=o("wks"),r=e.Symbol,l=c?r:r&&r.withoutSetter||s;n.exports=function(n){return f(u,n)||(u[n]=h&&f(r,n)?r[n]:l("Symbol."+n)),u[n]}},function(n,t,i){"use strict";var u=i(0),s=i(28).f,h=i(33),r=i(7),e=i(34),f=i(14),o=i(6),c=function(n){var t=function(t,i,r){if(this instanceof n){switch(arguments.length){case 0:return new n;case 1:return new n(t);case 2:return new n(t,i)}return new n(t,i,r)}return n.apply(this,arguments)};return t.prototype=n.prototype,t};n.exports=function(n,t){var it,a,y,i,l,p,w,k,d,v=n.target,g=n.global,rt=n.stat,ut=n.proto,b=g?u:rt?u[v]:(u[v]||{}).prototype,nt=g?r:r[v]||(r[v]={}),tt=nt.prototype;for(i in t)it=h(g?i:v+(rt?".":"#")+i,n.forced),a=!it&&b&&o(b,i),p=nt[i],a&&(n.noTargetGet?(d=s(b,i),w=d&&d.value):w=b[i]),l=a&&w?w:t[i],a&&typeof p==typeof l||(k=n.bind&&a?e(l,u):n.wrap&&a?c(l):ut&&"function"==typeof l?e(Function.call,l):l,(n.sham||l&&l.sham||p&&p.sham)&&f(k,"sham",!0),nt[i]=k,ut&&(y=v+"Prototype",!o(r,y)&&f(r,y,{}),r[y][i]=l,n.real&&tt&&!tt[i]&&f(tt,i,l)))}},function(n){var t={}.toString;n.exports=function(n){return t.call(n).slice(8,-1)}},function(n){n.exports=function(n){if(n==void 0)throw TypeError("Can't call method on "+n);return n}},function(n,t,i){var r=i(3),u=i(1),f=i(32);n.exports=!r&&!u(function(){return 7!=Object.defineProperty(f("div"),"a",{get:function(){return 7}}).a})},function(n,t,i){var r=i(3),u=i(15),f=i(4);n.exports=r?function(n,t,i){return u.f(n,t,f(1,i))}:function(n,t,i){return n[t]=i,n}},function(n,t,i){var f=i(3),e=i(13),r=i(36),o=i(5),u=Object.defineProperty;t.f=f?u:function(n,t,i){if(r(n),t=o(t,!0),r(i),e)try{return u(n,t,i)}catch(n){}if("get"in i||"set"in i)throw TypeError("Accessors not supported");return"value"in i&&(n[t]=i.value),n}},function(n,t,i){var r=i(1);n.exports=!!Object.getOwnPropertySymbols&&!r(function(){return!(Symbol()+"")})},function(n,t,i){var r,u,h=i(0),f=i(49),e=h.process,o=e&&e.versions,s=o&&o.v8;s?(r=s.split("."),u=r[0]+r[1]):f&&(r=f.match(/Edge\/(\d+)/),(!r||74<=r[1])&&(r=f.match(/Chrome\/(\d+)/),r&&(u=r[1])));n.exports=u&&+u},function(n,t,i){var r=i(7);n.exports=function(n){return r[n+"Prototype"]}},function(n,t,i){n.exports=i(23)},function(n,t,i){n.exports=i(51)},function(n,t,i){n.exports=i(22)},function(n,t,i){"use strict";i.r(t);var r=i(19),u=i.n(r),f=i(20),e=i.n(f);t["default"]={Render:function(n,t){function f(n,t){var r=document.getElementsByTagName("body")[0],i=document.createElement("script");i.src=n;t&&(i.onload=function(){return t(i)});r.appendChild(i)}function h(n){var i=document.getElementsByTagName("head")[0],t=document.createElement("link");t.rel="stylesheet";t.href=n;i.appendChild(t)}var i="",o,s,r;if(!t){if(r=e()(o=window.location.hostname.split(".")).call(o),2>r.length)throw"Invalid hostname ".concat(window.location.hostname);i=u()(s="".concat(r[1],".")).call(s,r[0]);t={userInfoUrl:"https://navigation-plugins-api.".concat(i,"/plugins-api/get-user-info"),baseUrl:"".concat(i),headerApiUrl:"https://navigation-plugins-api.".concat(i,"/plugins-api/get-header-resource")}}f(t.headerApiUrl,function(){h(headerNavigationPluginResource.info.StyleUrl);f(headerNavigationPluginResource.info.ScriptUrl,function(){headerNavigationPlugin.Render(n,t,headerNavigationPluginResource.menuConfiguration)})})}}},function(n,t,i){var r=i(24);n.exports=r},function(n,t,i){var u=i(25),r=Array.prototype;n.exports=function(n){var t=n.concat;return n===r||n instanceof Array&&t===r.concat?u:t}},function(n,t,i){i(26);var r=i(18);n.exports=r("Array").concat},function(n,t,i){"use strict";var o=i(10),s=i(1),h=i(8),c=i(2),l=i(37),a=i(38),r=i(40),v=i(41),y=i(48),p=i(9),w=i(17),u=p("isConcatSpreadable"),f=9007199254740991,e="Maximum allowed index exceeded",b=51<=w||!s(function(){var n=[];return n[u]=!1,n.concat()[0]!==n}),k=y("concat"),d=function(n){if(!c(n))return!1;var t=n[u];return t===void 0?h(n):!!t};o({target:"Array",proto:!0,forced:!b||!k},{concat:function(){for(var u,s,n,c=l(this),o=v(c,0),t=0,i=-1,h=arguments.length;i<h;i++)if(n=-1===i?c:arguments[i],d(n)){if(s=a(n.length),t+s>f)throw TypeError(e);for(u=0;u<s;u++,t++)u in n&&r(o,t,n[u])}else{if(t>=f)throw TypeError(e);r(o,t++,n)}return o.length=t,o}})},function(n){var t=function(){return this}();try{t=t||new Function("return this")()}catch(n){"object"==typeof window&&(t=window)}n.exports=t},function(n,t,i){var u=i(3),f=i(29),e=i(4),o=i(30),s=i(5),h=i(6),c=i(13),r=Object.getOwnPropertyDescriptor;t.f=u?r:function(n,t){if(n=o(n),t=s(t,!0),c)try{return r(n,t)}catch(n){}if(h(n,t))return e(!f.f.call(n,t),n[t])}},function(n,t){"use strict";var i={}.propertyIsEnumerable,r=Object.getOwnPropertyDescriptor,u=r&&!i.call({1:2},1);t.f=u?function(n){var t=r(this,n);return!!t&&t.enumerable}:i},function(n,t,i){var r=i(31),u=i(12);n.exports=function(n){return r(u(n))}},function(n,t,i){var r=i(1),u=i(11),f="".split;n.exports=r(function(){return!Object("z").propertyIsEnumerable(0)})?function(n){return"String"==u(n)?f.call(n,""):Object(n)}:Object},function(n,t,i){var f=i(0),u=i(2),r=f.document,e=u(r)&&u(r.createElement);n.exports=function(n){return e?r.createElement(n):{}}},function(n,t,i){var u=i(1),f=/#|\.prototype\./,r=function(n,t){var i=o[e(n)];return!(i!=h)||i!=s&&("function"==typeof t?u(t):!!t)},e=r.normalize=function(n){return(n+"").replace(f,".").toLowerCase()},o=r.data={},s=r.NATIVE="N",h=r.POLYFILL="P";n.exports=r},function(n,t,i){var r=i(35);n.exports=function(n,t,i){return(r(n),void 0===t)?n:0===i?function(){return n.call(t)}:1===i?function(i){return n.call(t,i)}:2===i?function(i,r){return n.call(t,i,r)}:3===i?function(i,r,u){return n.call(t,i,r,u)}:function(){return n.apply(t,arguments)}}},function(n){n.exports=function(n){if("function"!=typeof n)throw TypeError(n+" is not a function");return n}},function(n,t,i){var r=i(2);n.exports=function(n){if(!r(n))throw TypeError(n+" is not an object");return n}},function(n,t,i){var r=i(12);n.exports=function(n){return Object(r(n))}},function(n,t,i){var r=i(39),u=Math.min;n.exports=function(n){return 0<n?u(r(n),9007199254740991):0}},function(n){var t=Math.ceil,i=Math.floor;n.exports=function(n){return isNaN(n=+n)?0:(0<n?i:t)(n)}},function(n,t,i){"use strict";var r=i(5),u=i(15),f=i(4);n.exports=function(n,t,i){var e=r(t);e in n?u.f(n,e,f(0,i)):n[e]=i}},function(n,t,i){var u=i(2),r=i(8),f=i(9),e=f("species");n.exports=function(n,t){var i;return r(n)&&(i=n.constructor,"function"==typeof i&&(i===Array||r(i.prototype))?i=void 0:u(i)&&(i=i[e],null===i&&(i=void 0))),new(void 0===i?Array:i)(0===t?0:t)}},function(n,t,i){var u=i(43),r=i(44);(n.exports=function(n,t){return r[n]||(r[n]=t===void 0?{}:t)})("versions",[]).push({version:"3.6.2",mode:u?"pure":"global",copyright:"© 2020 Denis Pushkarev (zloirock.ru)"})},function(n){n.exports=!0},function(n,t,i){var u=i(0),f=i(45),r="__core-js_shared__",e=u[r]||f(r,{});n.exports=e},function(n,t,i){var r=i(0),u=i(14);n.exports=function(n,t){try{u(r,n,t)}catch(i){r[n]=t}return t}},function(n){var t=0,i=Math.random();n.exports=function(n){return"Symbol("+((n===void 0?"":n)+"")+")_"+(++t+i).toString(36)}},function(n,t,i){var r=i(16);n.exports=r&&!Symbol.sham&&"symbol"==typeof Symbol.iterator},function(n,t,i){var r=i(1),u=i(9),f=i(17),e=u("species");n.exports=function(n){return 51<=f||!r(function(){var t=[],i=t.constructor={};return i[e]=function(){return{foo:1}},1!==t[n](Boolean).foo})}},function(n,t,i){var r=i(50);n.exports=r("navigator","userAgent")||""},function(n,t,i){var r=i(7),u=i(0),f=function(n){if("function"==typeof n)return n};n.exports=function(n,t){return 2>arguments.length?f(r[n])||f(u[n]):r[n]&&r[n][t]||u[n]&&u[n][t]}},function(n,t,i){var r=i(52);n.exports=r},function(n,t,i){var u=i(53),r=Array.prototype;n.exports=function(n){var t=n.reverse;return n===r||n instanceof Array&&t===r.reverse?u:t}},function(n,t,i){i(54);var r=i(18);n.exports=r("Array").reverse},function(n,t,i){"use strict";var u=i(10),f=i(8),e=[].reverse,r=[1,2];u({target:"Array",proto:!0,forced:r+""==r.reverse()+""},{reverse:function(){return f(this)&&(this.length=this.length),e.call(this)}})}])["default"]})