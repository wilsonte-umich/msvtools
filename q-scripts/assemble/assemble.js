function getUrlVars(){
    var vars = [], hash;
    var hashes = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');
    for(var i = 0; i < hashes.length; i++)
    {
        hash = hashes[i].split('=');
        vars.push(hash[0]);
        vars[hash[0]] = hash[1];
    }
    return vars;
}
function sortList(id) {
  var list, i, switching, b, shouldSwitch;
  list = document.getElementById(id);
  switching = true;
  while (switching) {
    switching = false;
    b = list.getElementsByTagName("li");
    for (i = 0; i < (b.length - 1); i++) {
      shouldSwitch = false;
      if (b[i].innerHTML.toLowerCase() > b[i + 1].innerHTML.toLowerCase()) {
        shouldSwitch = true;
        break;
      }
    }
    if (shouldSwitch) {
      b[i].parentNode.insertBefore(b[i + 1], b[i]);
      switching = true;
    }
  }
}
//---------------------------------------------------------
// jquery
//---------------------------------------------------------
$(function() { // initialize the main page
    pge.init();
});
//---------------------------------------------------------
// ajax communication with server
//---------------------------------------------------------
var pge = { // generic functions for ajax server calls
    sel: 'selected-sample',
    crr: undefined,
    init: function(){ // general function for debugging
        var qry = getUrlVars();
        
        console.log(qry);
        
        if(qry.blinded == 0){
            $.each($("#sample-list li"), function(i, li){
                $(li).html($(li).data("sample"));
            });
            sortList('sample-list');
        }
        $("#sample-list li").on("click", function(){
            sel.do($(this));
        });
        $("#selector-buttons button").on("click", function(){
            if (pge.crr == undefined){
                sel.do($("#sample-list li")[1]);
            } else if($(this).data("inc") > 0) {
                sel.do(pge.crr.next());
            } else {
                sel.do(pge.crr.prev());
            }   
        });  
    }
};
var sel = {
    do: function(li){
        if(li.length > 0){
            var smp = li.data("sample"),
                pfx = "samples/plots/" + smp + "/" + smp;
            $("#correlation-plot").attr('src', pfx + '.CN_OUT.CORRELATION.jpg');
            $("#cn-change-plot")  .attr('src', pfx + '.CN_OUT.CN_IN.CN.jpg');
            $("#lrr-plot")        .attr('src', pfx + '.CN_OUT.CN_2_LRR.jpg');
            $("#sample-list li").removeClass(pge.sel);
            li.addClass(pge.sel);
            $("#current-sample").html(smp);
            pge.crr = li;               
        }
    }
}


