library(networkD3)
datf = 
  read.csv("file_with_edge_weight_and_identity_data.csv", 
           header = TRUE)
datf$Edge_weight = datf$Edge_weight*100
head(datf, n = 3)
names = 
  read.csv("file_with_identity_and_class_data.csv",
           header = TRUE)
head(names, n = 3)
datf$Edge_type <- as.factor(datf$Edge_type)
names$Ancestor_type <- as.factor(names$Ancestor_type)
my_color <- 
  'd3.scaleOrdinal() .range([ 
"#4F4F4F","#696969","#909090","#cbcbcb",
"#8DD3C7","#D9D9D9","#4169E1","#8cb8d7","#BEBADA","#FDB462","#FB8072",
"#FCCDE5","#FFFFB3"])'
sankey = sankeyNetwork(Links = datf, Nodes = names, Source = "Number_From",
                       Target = "Number_To", Value = "Edge_weight", 
                       NodeID = "cell_states", units = "Edge weight",
                       colourScale = my_color, 
                       LinkGroup="Edge_type", NodeGroup="Ancestor_type",
                       fontSize = 36, nodeWidth = 40, iterations = 0,
                       sinksRight = FALSE,
                       height = 3000, width = 3000)
sankey
sankey = sankeyNetwork(Links = datf, Nodes = names, Source = "Number_From",
                       Target = "Number_To", Value = "Edge_weight", 
                       NodeID = "cell_states", units = "Edge weight",
                       colourScale = my_color, 
                       LinkGroup="Edge_type", NodeGroup="Ancestor_type",
                       fontSize = 0, nodeWidth = 50, iterations = 0,
                       sinksRight = FALSE,
                       height = 6000, width = 6000, margin = list(bottom = 150))
sankey
library(htmlwidgets)
onRender(sankey,
         '
  function(el, x) {
    var sankey = this.sankey;
    var path = sankey.link();
    var nodes = d3.selectAll(".node");
    var link = d3.selectAll(".link")
    var width = el.getBoundingClientRect().width - 40;
    var height = el.getBoundingClientRect().height - 40;

    window.dragmove = function(d) {
      d3.select(this).attr("transform", 
        "translate(" + (
           d.x = Math.max(0, Math.min(width - d.dx, d3.event.x))
            ) + "," + (
            d.y = Math.max(0, Math.min(height - d.dy, d3.event.y))
          ) + ")");
      sankey.relayout();
      link.attr("d", path);
    };

    nodes.call(d3.drag()
      .subject(function(d) { return d; })
      .on("start", function() { this.parentNode.appendChild(this); })
      .on("drag", dragmove));
  }
  '
)
setwd("your_directory")
saveNetwork(sankey, "sankey_all.html")
#install.packages("webshot2")
library(webshot2)
chromote::set_chrome_args("--disable-crash-reporter")
webshot2::webshot("sankey_all_labelled.html", 
                  "sankey_all_labelled.png", vwidth = 3000, vheight = 3000)
