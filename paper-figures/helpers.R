CreatePair = function(x,y, collapse = "&") {paste(sort(c(unique(x), unique(y))), collapse = collapse)}
CreatePair = Vectorize(CreatePair)

CheckAAcode = function(cds) {
  aasequence = as.character(cds)
  unknown_chars_in_seq = gsub(aalettersstr, "", aasequence)
  return(nchar(unknown_chars_in_seq))
}

# antidefense hits
GetAccessionShort = function(accession) {
  gsub(" ", "", accession)
  strsplit(accession, split = "\\.")[[1]][1]
}


Parse.Hhr.Raw.Result = function(hhblits.table) {
  hhblits.table = hhblits.table[,qstart := as.numeric(qstart)]
  hhblits.table = hhblits.table[,qend := as.numeric(qend)]
  hhblits.table = hhblits.table[,sstart := as.numeric(sstart)]
  hhblits.table = hhblits.table[,send := as.numeric(send)]
  hhblits.table = hhblits.table[,prob := as.numeric(prob)]
  hhblits.table = hhblits.table[,eval := as.numeric(eval)]
  hhblits.table = hhblits.table[,qcov := round((qend - qstart + 1)/qlength,d=4)]
  hhblits.table = hhblits.table[,scov := round((send - sstart + 1)/slength,d=4)]
  hhblits.table = hhblits.table[,annotation := tolower(annot)]
  hhblits.table = hhblits.table[,hit.length := qend - qstart + 1]
  return(hhblits.table)
}


Get.uid.From.ECOD.subject.id = function(subject.id) {
  uid = strsplit(subject.id, "_")[[1]][2]
  return(uid)
}


Get.ecod.domain.id.From.ECOD.subject.id = function(subject.id) {
  ecod.domain.id = strsplit(subject.id, "_")[[1]][3]
  return(ecod.domain.id)
}

Annotate_Graph = function(data.for.network, 
                               edge.colors = NULL,
                               edge.linetypes = NULL,
                               node.color.map,
                               node.size.map = NULL,
                               node.labels = NULL,
                               line.color.variable,
                               vertex.size = 1, 
                               remove.loops = TRUE) {
  if (is.null(node.size.map)) {
    node.size.map = data.frame(node = unique(c(data.for.network$from, data.for.network$to))) %>% mutate(size =vertex.size)
  }
  if (is.null(node.labels)) {
    node.labels = data.frame(node = unique(c(data.for.network$from, data.for.network$to))) %>% mutate(label = node)
  }
  label.color.map = node.color.map
  data.for.network$edge.var =  data.for.network %>% pull(line.color.variable)
  net_hh = igraph::graph_from_data_frame(data.for.network)# %>% filter(min.cov == 0.3))
  df_edges <- igraph::as_data_frame(net_hh, what = "edges") %>% 
    select(from, to, edge.var)
  
  net_hh = igraph::simplify(net_hh,
                              remove.multiple = TRUE,
                              remove.loops = remove.loops)
    
    df_edges_simpliied <- igraph::as_data_frame(net_hh, what = "edges")
    df_edges <- df_edges_simpliied %>% 
      left_join(df_edges %>% select(from, to, edge.var)) 

  #left_join(data.for.network %>% distinct(from, to, edge.var))
  new_graph <- igraph::graph_from_data_frame(d = df_edges, vertices = igraph::as_data_frame(net_hh, what = "vertices"))
  nodes =  igraph::V(new_graph)
  net_annotated = new_graph
  for (col in unique(node.color.map$color)) {
    seqs = node.color.map %>% filter(color == col) %>% distinct(node) %>% pull(node)
    vertex = nodes[names(nodes)  %in% seqs]
    net_annotated = net_annotated %>%
      igraph::set_vertex_attr("color", index = vertex, col)
  }
  for (this.size in unique(node.size.map$size)) {
    seqs = node.size.map %>% filter(size == this.size) %>% distinct(node) %>% pull(node)
    vertex = nodes[names(nodes)  %in% seqs]
    net_annotated = net_annotated %>%
      igraph::set_vertex_attr("size", index = vertex, this.size)
  }
  for (this.col in unique(label.color.map$color)) {
    seqs = label.color.map %>% filter(color ==this.col) %>% distinct(node) %>% pull(node)
    vertex = nodes[names(nodes)  %in% seqs]
    net_annotated = net_annotated %>%
      igraph::set_vertex_attr("label.color", index = vertex, this.col)
  }
  #net_annotated = net_annotated %>% igraph::simplify()
  if (is.null(edge.colors)) {
    igraph::E(net_annotated)$edge.color = "gray"
  } else {
    for (i in 1:nrow(edge.colors)) {
      igraph::E(net_annotated)$edge.color[igraph::E(net_annotated)$edge.var == edge.colors$value[i]] <- edge.colors$color[i]
    }
  }
  
  if (is.null(edge.linetypes)) {
    igraph::E(net_annotated)$lty = "solid"
  } else {
    for (i in 1:nrow(edge.linetypes)) {
      igraph::E(net_annotated)$lty[igraph::E(net_annotated)$edge.var == edge.linetypes$value[i]] <- edge.linetypes$lty[i]
    }
  }
  #net_annotated = net_annotated %>% igraph::simplify()
    for (this.label in unique(node.labels$label)) {
      seqs = node.labels %>% filter(label == this.label) %>% distinct(node) %>% pull(node)
      vertex = nodes[names(nodes)  %in% seqs]
      net_annotated = net_annotated %>%
        igraph::set_vertex_attr("label.name", index = vertex, this.label)
    }
  return(net_annotated)
}


#' Visualie network using plotly package.
#' We need to first convert our network to a data.frame using ggnetwork package
#' @net_annotated network with annotation attribute
#' @colors a mapping of annotations and colors
VisualiseGGnetwork = function(net_annotated, vertex.size = 0.5, cex = 0.1, 
                              label_size = 1, 
                              layout = igraph::nicely(), 
                              main = "", 
                              nudge_y = 0.025, 
                              alpha = 1,
                              remove.loops = T,
                              simplify.net = TRUE) {
  
  if (simplify.net == TRUE) {
  data = ggnetwork(x = igraph::simplify(net_annotated, remove.multiple = T, remove.loops = remove.loops), 
                   layout = layout,
                   scale = TRUE, 
                   stringsAsFactors = FALSE, 
                   arrow.gap = 0)
  }
  
  
  if (!("size" %in% names(data))) {
    data$size = vertex.size
  }
  if (!("color" %in% names(data))) {
    data$color = "gray"
  }
  if(!("label.name" %in% names(data))) {
    data$label.name = data$name
  }
  
  g=ggplot(data) +
    geom_edges(aes(x = x, y = y, xend=xend, yend=yend), color = "grey50", size = cex) +
    geom_nodes(aes(x = x, y = y, color = color, label1 = label.name, size = size), alpha = alpha) +
    geom_text(aes(x = x, y = y, color = color, label = label.name), size = label_size, nudge_y = nudge_y) +
    guides(size = "none") +
    scale_color_identity() +
    xlim(c(-0.05, 1.05)) +
    theme_blank() +
    ggtitle(main)
  
  return(g)
}


Get.Ecod.Annotation.Mapping = function(distinct_domains, ecod.domain.descriptions) {

  unclassified.f.names = c("F_UNCLASSIFIED", "UNK_F_TYPE")
  ecod.domain.descriptions.full.names = ecod.domain.descriptions %>%
    rowwise() %>%
    tidyr::separate(col = "f_id", into = c('x_ind', 'h_ind', 't_ind', 'f_ind'), sep = "\\.", remove = FALSE) %>%
    # we will have no f_ind in some cases
    mutate(f_unclassified = f_name %in% unclassified.f.names) 
  
  # WE COULD ADD F_NAME where it is known in cases when one number is missing but in our data there are actually no cases of known f_name for missing f number
  #idx.uknown.f.id.with.known.f.name = which(is.na(ecod.domain.descriptions.full.names$f_ind) & !ecod.domain.descriptions.full.names$f_unclassified)
  #ecod.domain.descriptions.full.names$f_ind[idx.uknown.f.id.with.known.f.name] = ecod.domain.descriptions.full.names$f_name[idx.uknown.f.id.with.known.f.name]
  
  ecod.domain.descriptions.full.names = ecod.domain.descriptions.full.names  %>%
    rename(f_id_orig = f_id) %>%
    mutate(x_id = x_ind) %>%
    tidyr::unite(col = 'h_id', x_ind, h_ind, sep = ".", remove = FALSE)  %>%
    tidyr::unite(col = 't_id', x_ind, h_ind, t_ind, sep = ".", remove = FALSE) %>%
    tidyr::unite(col = 'f_id', x_ind, h_ind, t_ind, f_ind, sep = ".", remove = FALSE) 

  ecod.annotation.map = data.frame(domain = distinct_domains, stringsAsFactors = FALSE)
  ecod.annotation.map$X.uid =  sapply(1:nrow(ecod.annotation.map),function(row.index){Get.uid.From.ECOD.subject.id(ecod.annotation.map$domain[row.index])})
  ecod.annotation.map$ecod_domain_id = sapply(1:nrow(ecod.annotation.map),function(row.index){Get.ecod.domain.id.From.ECOD.subject.id(ecod.annotation.map$domain[row.index])})
  ecod.annotation.map = ecod.annotation.map %>%
    left_join(ecod.domain.descriptions.full.names, by = c("X.uid", "ecod_domain_id"))
  return(ecod.annotation.map)
}



Get.Ecod.Id.To.Name.Map = function(ecod.annotation.map) {
  # now find the names for each id
  ecod_id_to_names_map = ecod.annotation.map %>% 
    group_by(x_id) %>% 
    summarise(name = paste(sort(unique(x_name)), collapse = " | ")) %>% 
    select(id = x_id, name) %>%
    mutate(level = "X")  %>%
    rbind(ecod.annotation.map %>% 
            group_by(h_id) %>% 
            summarise(name = paste(sort(unique(h_name)), collapse = " | ")) %>% 
            select(id = h_id, name) %>%
            mutate(level = "H")) %>%
    rbind(ecod.annotation.map %>% 
            group_by(t_id) %>% 
            summarise(name = paste(sort(unique(t_name)), collapse = " | ")) %>% 
            select(id = t_id, name) %>%
            mutate(level = "T") ) %>%
    rbind(ecod.annotation.map %>% 
            group_by(f_id) %>% 
            summarise(name = paste(sort(unique(f_name)), collapse = " | ")) %>% 
            select(id = f_id, name) %>%
            mutate(level = "F") )
  
  
  no_h_name_map = ecod.annotation.map %>% 
    filter(h_name == "NO_H_NAME") %>% 
    group_by(h_id) %>% 
    summarise(h_name = paste(unique(t_name), collapse = " | ", sep = "")) %>%
    select(id = h_id, name = h_name) %>%
    mutate(level = "H")
  
  ecod_id_to_names_map = rbind(ecod_id_to_names_map %>% filter(name != "NO_H_NAME"),
                               no_h_name_map)
  return(ecod_id_to_names_map)
}



Get.Domain.Positions.Data = function(annotated.ecod.domains) {
  annotated.ecod.domains = annotated.ecod.domains %>% 
    select(qname, annotation, qlength, qstart, qend, domain, family)
  tile.data = data.frame(qname = character(0), 
                         pos.start = integer(0), 
                         pos.end = integer(0), 
                         domain = character(0))
  if (nrow(annotated.ecod.domains) == 0) {
    return(tile.data)
  } else {
    for (this.qname in unique(annotated.ecod.domains$qname)) {
      dat = annotated.ecod.domains %>% filter(qname == this.qname)
      qlength = dat %>% pull(qlength) %>% unique()
      all.borders = sort(unique(c(1, qlength, dat$qstart, dat$qend)))
      for (i in seq(2,length(all.borders))) {
        left.boarder = all.borders[i-1]
        right.boarder = all.borders[i]
        domain.this.part = dat %>% filter(qstart <= left.boarder & qend >= right.boarder) %>% pull(domain) %>% unique() %>% sort()
        this.qname.this.part.data = data.frame(qname = this.qname, 
                                               pos.start = if_else(left.boarder %in% c(1, dat$qstart), left.boarder, left.boarder +1), 
                                               pos.end = if_else(right.boarder %in% c(dat$qend, qlength), right.boarder, right.boarder-1), 
                                               domain =  paste(domain.this.part, sep = " & ", collapse = " & "),
                                               multi_domain = length(domain.this.part)>1)
        tile.data = rbind(tile.data, this.qname.this.part.data)
      }
    }
    tile.data = tile.data %>%
      mutate(
        domain = if_else(domain == "", "undetected", domain),
        domain = if_else(multi_domain, "multiple domains", domain)) %>%
      rowwise() %>%
      mutate(row = which(unique(annotated.ecod.domains$qname) == qname) -0.5) %>%
      mutate(row.plus = row + 0.5) %>%
      mutate(qname = factor(qname, levels = unique(annotated.ecod.domains$qname)))
    return(tile.data)
  }
}


Show.Domain.Position.Within.Data = function(tile.data, colormap, my.theme = theme_bw(), 
                                            color.frames = FALSE, legend.position = "bottom",
                                            interactive.plot = TRUE) {
  if (nrow(tile.data) == 0) {
    gg_rect = ggplot() + geom_blank()
  } else {
    tile.data = tile.data %>%
      mutate(tooltip = paste0("Domain annotation:", domain_name, "\n", "seq: ", qname, "\n", "domain.id: ", domain)) %>%
      arrange(row)
    
    tile.data$qname = factor(tile.data$qname, levels = unique(tile.data$qname))
    tile.data$family = factor(tile.data$family, levels = unique(tile.data$family))
    
    if (color.frames) {
      gg_rect = ggplot(tile.data, aes(color = x_name))
    } else {
      gg_rect0 = ggplot(tile.data, aes(color = domain))  +
        scale_color_manual(values = colormap[which(names(colormap) %in% unique(tile.data$domain))], name = "Domain")
    }
    if (interactive.plot) {
      gg_rect1 = gg_rect0 +   
        geom_blank(aes(x=pos.start, y = qname)) +
        geom_rect_interactive(mapping = aes(xmin = pos.start, xmax = pos.end,
                                            ymin = row, ymax = row.plus, fill = domain, 
                                            tooltip = tooltip), alpha=1)       
    } else {
      gg_rect1 = gg_rect0 +   
        geom_blank(aes(x=pos.start, y = qname)) +
        geom_rect(mapping = aes(xmin = pos.start, xmax = pos.end,
                                            ymin = row, ymax = row.plus, fill = domain), alpha=1) 
    }
  gg_rect = gg_rect1 +
      scale_x_continuous(name="position", limits = c(0, max(tile.data$pos.end)),  expand = c(0, 0)) +
      scale_y_discrete(name = "", expand = c(0, 0)) +
      # leave in the egend only the colors that we can actually see
      scale_fill_manual(values = colormap[which(names(colormap) %in% unique(tile.data$domain))], 
                        name = "Domain",
                        guide = guide_legend(ncol = 2)) + 
      facet_grid(. ~annotation) + 
      my.theme +
      theme(legend.position = legend.position) #+
      #guides(
      #  color = guide_legend(show = FALSE),
      #  fill = guide_legend(nrow = 5,byrow=TRUE)
      #)
  }
  #guides(fill = "none")
  return(gg_rect)
}

# Set ggplot theme
Get.Theme= function(text_size = 8, angle.x = 45, angle.y=0, legend.position = "left") {
  my_theme = theme(
    panel.grid.major = element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size=text_size),
    axis.text.x = element_text(angle=angle.x, hjust = 1,vjust=1, size = text_size),
    axis.text.y = element_text(angle=angle.y, hjust = 1,vjust=1, size = text_size),
    legend.position = legend.position,
    panel.border = element_rect(linetype = "solid", fill = NA, color = "black"))
  return(my_theme)
} 


Get.significant.phrog.hits = function(hhr.phrogs, min.prob.threshold, min.scov.threshold, min.qcov.threshold, max.eval.threshold) {
  significant.hhr.phrogs = hhr.phrogs %>%
    filter(prob >= min.prob.threshold,
           scov >= min.scov.threshold,
           qcov >= min.qcov.threshold,
           eval <= max.eval.threshold)
  return(significant.hhr.phrogs)
}



GetPhrogAnnotationsData = function(hhr.phrogs.annotated, annotation.indexes, min.cov, min.prob, min.eval) {
  annotated.proteins.including.multi.annot = Get.significant.phrog.hits(hhr.phrogs.annotated, 
                                                                        min.prob.threshold = min.prob, 
                                                                        min.qcov.threshold = min.cov, 
                                                                        min.scov.threshold = min.cov, 
                                                                        max.eval.threshold = min.eval) %>%
    distinct(qname, annotation, category, include) %>%
    left_join(annotation.indexes, by = c("annotation", "category")) %>% 
    select(qname, annotation.index, annotation, category, include) %>%
    group_by(qname) %>%
    # number of annotations including the excluded ones
    mutate(num.annots = n_distinct(annotation),
           num.antidefense.annots = n_distinct(annotation[category == "antidefense"])) %>%
    ungroup() 
  return(annotated.proteins.including.multi.annot)
}


GetUniquelyAnnotatedProteins = function(annotated.proteins.including.multi.annot, GENERAL.ANTIDEFENSE.ANNOTATION.INDEX, unspecific.annotation.indices) {
 
  antidefense.proteins = annotated.proteins.including.multi.annot %>%
    filter(num.antidefense.annots > 0) %>%
    # only keep the antidefense annotations as they are more explicit
    filter(category == "antidefense") %>%
    group_by(qname, family) %>%
    summarise(annotation.index = if_else(num.antidefense.annots > 1,  as.integer(GENERAL.ANTIDEFENSE.ANNOTATION.INDEX), annotation.index),
              include = if_else(num.antidefense.annots > 1,  TRUE, include),
              num.annots = 1) %>%
    ungroup() %>%
    select(qname, family, annotation.index, include, num.annots)
  
  
   unspecific.proteins.with.multiple.annotations = annotated.proteins.including.multi.annot %>%
      # the ones with antidefense annot were dealed with above
      inner_join(unspecific.annotation.indices) %>%
      filter(num.antidefense.annots == 0 & num.annots > 1) %>%
      distinct(qname) %>%
      left_join(annotated.proteins.including.multi.annot) %>%
      # get rid of unspecific annotations
      anti_join(unspecific.annotation.indices) %>%
      group_by(qname) %>%
      # now we count the number of specific annotations 
      mutate(num.annots = n_distinct(annotation.index)) %>%
      ungroup() %>%
      select(qname, family, annotation.index, include, num.annots)
  

  annotated.proteins = 
    annotated.proteins.including.multi.annot %>% 
    anti_join(antidefense.proteins %>% select(qname), by = "qname") %>%
    anti_join(unspecific.proteins.with.multiple.annotations %>% select(qname), by = "qname")%>%
    select(qname, family, annotation.index, include, num.annots) %>%
    rbind(antidefense.proteins) %>%
    rbind(unspecific.proteins.with.multiple.annotations) %>%
    filter(num.annots == 1) %>%
    select(qname, annotation.index, family, include)

  return(annotated.proteins)
}



Calculate.Annotation.Tradeoff = function(hhr.phrogs.annotated, max.eval.range, min.cov.range, min.prob.range, num.proteins) {
  for (this.eval in max.eval.range)  {
    for (this.cov in min.cov.range)  {
      for (this.prob in min.prob.range) {
        annotation.uncertainty.data = Get.significant.phrog.hits(hhr.phrogs.annotated, this.prob, this.cov, this.cov, this.eval)  %>%
          group_by(qname) %>%
          summarise(num.annotations = n_distinct(annotation),
                    num.categories = n_distinct(category))  %>%
          summarise(num.prot = num.proteins,
                    num.annotated = n(),
                    perc.annotated = 100*num.annotated/num.proteins,
                    perc.unique.annotation = 100*sum(num.annotations == 1)/num.annotated) %>%
          #perc.unique.category = 100*sum(num.categories == 1)/num.annotated) 
          mutate(annotation.coverage = this.cov, minimum.probability = this.prob, maximum.eval = this.eval) %>%
          select(annotation.coverage,minimum.probability, maximum.eval, num.prot, num.annotated, perc.annotated, perc.unique.annotation) %>%
          rbind(annotation.uncertainty.data)
      }
    }
  }
  return(annotation.uncertainty.data)
}

VisualiseIgraphnetwork = function(data.for.network, 
                               edge.colors,
                               node.color.map,
                               node.size.map = NULL,
                               node.labels,
                               edge.linetypes = NULL,
                               line.color.variable,
                               line.color.legend.name,
                               font.size = 1,
                               vertex.size = 1,
                               add.legend = TRUE,
                               main = "",
                               vertex.label.dist = 0.6,
                               include.labels = TRUE, 
                               remove.loops = FALSE,
                               alpha = 1,
                               label.font = 1, 
                               legend.title = NULL,
                               node.label.color = NULL) {  
  
    
  net_annotated =  Annotate_Graph(data.for.network, 
                                       edge.colors = edge.colors,
                                       edge.linetypes = edge.linetypes,   
                                       node.color.map = node.color.map,
                                       node.size.map = node.size.map,
                                       node.labels = node.labels,
                                       line.color.variable = line.color.variable,
                                       vertex.size = vertex.size,
                                       remove.loops = remove.loops)
  if (include.labels) {
    vertex.label = igraph::vertex_attr(net_annotated, "name")
  } else {
    vertex.label = NA 
  }
  
  if (is.null(igraph::E(net_annotated)$weight)) {
    igraph::E(net_annotated)$weight = 0.1
  }
  if (is.null(node.label.color)) {
    node.label.color = igraph::V(net_annotated)$color
  }
  
  p.plot=igraph::plot.igraph(
    net_annotated,
    labels = igraph::vertex_attr(net_annotated, "name"),
    vertex.label = igraph::V(net_annotated)$label.name,
    #igraph::vertex_attr(net_annotated, "label"),
    vertex.label.cex = font.size,
    vertex.label.font = 1,
    vertex.label.dist = vertex.label.dist,
    vertex.label.color = node.label.color,
    vertex.frame.color = adjustcolor(igraph::V(net_annotated)$color,alpha.f = alpha),
    vertex.color = adjustcolor(igraph::V(net_annotated)$color,alpha.f = alpha),
    vertex.edge.color =  igraph::V(net_annotated)$color,
    vertex.size =  igraph::V(net_annotated)$size, 
    arrow.size = 0.2, 
    edge.arrow.size = 0,
    main = main,
    layout =  igraph::layout_nicely,
    edge.width = 10*igraph::E(net_annotated)$weight,
    edge.color = igraph::E(net_annotated)$edge.color,
    edge.lty = igraph::E(net_annotated)$lty,
    cex.main=5,
    label.font = label.font
    #edge.curved = igraph::curve_multiple(net_annotated, .2)
  )
  if (!is.null(edge.colors) & add.legend) {
    edges = inner_join(edge.linetypes, edge.colors, by = "value")
    legend(1,1,
           title = legend.title, 
           col = edges$color,
           lty = edges$lty,
           #legend = paste0(line.color.legend.name," = ", edge.colors$value),
           legend = edges$value,
           text.font = 1,
           cex = font.size)
  }
  return(p.plot)
}

Get.Pairs.That.Share.T = function(ecods.across.fold.types.data) {
  ecods.across.fold.types.data.stripped = ecods.across.fold.types.data %>%
    select(t_id, fold.type, all.x, num.x)
  
  fold.type.pairs.that.share.t = ecods.across.fold.types.data.stripped %>%
    inner_join(ecods.across.fold.types.data.stripped, by = c("t_id"), relationship = "many-to-many") %>%
    filter(fold.type.x != fold.type.y) %>%
    rename('shared_t_id' = 't_id')
  return(fold.type.pairs.that.share.t)
}

Get.Pairs.That.Share.T.And.Disshare.X = function(ecods.across.fold.types.data) {
  fold.type.pairs.that.share.t = Get.Pairs.That.Share.T(ecods.across.fold.types.data) 
  fold.type.pairs.that.share.t.and.disshare.x = fold.type.pairs.that.share.t %>%
    filter(num.x.x > 1 & num.x.y > 1) %>%
    # now we have fold type pairs that share t
    # now to each pair add all x and t seen in this fold type
    # so we will have more rows now: multiple rows per each fold type
    left_join(ecods.across.fold.types %>% distinct(fold.type, x_id, all.x, annotation.index) %>% select(fold.type.x = fold.type, x.id.x = x_id, annotation.index.x = annotation.index, all.x.within.fold.type.x =all.x)) %>%
    left_join(ecods.across.fold.types %>% distinct(fold.type, x_id, all.x, annotation.index) %>% select(fold.type.y = fold.type, x.id.y = x_id, annotation.index.y = annotation.index, all.x.within.fold.type.y =all.x)) %>%
    group_by(fold.type.x, fold.type.y, shared_t_id) %>%
    mutate(all.x.within.both.fold.types = paste(sort(unique(c(x.id.x, x.id.y))), collapse = " & ")) %>%
    ungroup() %>%
    filter(all.x.within.both.fold.types != all.x.within.fold.type.x & all.x.within.both.fold.types != all.x.within.fold.type.y) %>%
    distinct(shared_t_id, fold.type.x, fold.type.y, shared_t_id, annotation.index.x, annotation.index.y)
  return(fold.type.pairs.that.share.t.and.disshare.x)
}

Calculate.Num.Mosaic.Fold.Types = function(ecods.across.fold.types.data) {
  fold.type.pairs.that.share.t.and.disshare.x = Get.Pairs.That.Share.T.And.Disshare.X(ecods.across.fold.types.data) 
   mosaic.pairs.num = fold.type.pairs.that.share.t.and.disshare.x %>%
    distinct(fold.type.x, fold.type.y, annotation.index.x, annotation.index.y) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(fold.type.x, fold.type.y)), collapse = " AND ", sep = "")) %>%
    group_by(annotation.index.x, annotation.index.y)  %>%
   summarise(num.mosaic.pairs = n_distinct(pair)) %>%
   ungroup()
return(mosaic.pairs.num)
}




Calculate.Mosaicism.Per.Fold.Type = function(ecods.across.fold.types.data) {
  num.pairs.of.fold.types = ecods.across.fold.types.data %>% 
    group_by(annotation.index) %>%
    summarise(num.fold.types = n_distinct(fold.type),
              num.fold.types.with.multi.x = n_distinct(fold.type[num.x > 1])) %>%
    ungroup() %>%
    mutate(num.theoretical.fold.type.pairs = num.fold.types*(num.fold.types-1)/2,
           num.theoretical.fold.type.pairs.with.multi.x = num.fold.types.with.multi.x*(num.fold.types.with.multi.x-1)/2)
  
  # now look only within category
  fold.type.pairs.that.share.t = Get.Pairs.That.Share.T(ecods.across.fold.types.data) %>%
    inner_join(ecods.across.fold.types %>% distinct(fold.type, annotation.index) %>% select(fold.type.x = fold.type, annotation.index)) %>%
    # join by annotation.index so that we only look within groups
    inner_join(ecods.across.fold.types %>% distinct(fold.type, annotation.index) %>% select(fold.type.y = fold.type, annotation.index)) 
  
  fold.type.pairs.that.share.t.and.disshare.x = Get.Pairs.That.Share.T.And.Disshare.X(ecods.across.fold.types.data) %>%
    inner_join(ecods.across.fold.types %>% distinct(fold.type, annotation.index) %>% select(fold.type.x = fold.type, annotation.index)) %>%
    # join by annotation.index so that we only look within groups
    inner_join(ecods.across.fold.types %>% distinct(fold.type, annotation.index) %>% select(fold.type.y = fold.type, annotation.index)) 
  
    
  all.fold.type.pairs.sharing.t = fold.type.pairs.that.share.t %>%
    distinct(fold.type.x, fold.type.y, annotation.index) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(fold.type.x, fold.type.y)), collapse = " AND ", sep = "")) %>% 
    group_by(annotation.index)  %>%
    summarise(num.pairs.sharing.t = n_distinct(pair)) %>%
    ungroup()
  
 
 mosaic.pairs.num = fold.type.pairs.that.share.t.and.disshare.x %>%
    distinct(fold.type.x, fold.type.y, annotation.index) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(fold.type.x, fold.type.y)), collapse = " AND ", sep = "")) %>% 
    group_by(annotation.index)  %>%
    summarise(num.mosaic.pairs = n_distinct(pair)) %>%
    ungroup()
    
  mosaicism.per.fold.type   = ecods.across.fold.types.data %>% 
    distinct(annotation.index) %>%
    left_join(mosaic.pairs.num) %>%
    left_join(all.fold.type.pairs.sharing.t) %>%
    left_join(num.pairs.of.fold.types) %>%
    mutate(prob.mosaicism = num.mosaic.pairs/num.pairs.sharing.t,
           prop.pairs.mosaoc.by.domain = num.mosaic.pairs/num.theoretical.fold.type.pairs)
  
  return(mosaicism.per.fold.type)
}


Get.Num.Mosaic.Family.Pairs = function(all.family.pairs.sharing.t.and.disshare.x, all.family.pairs.sharing.t) {
  mosaicism.per.family = all.family.pairs.sharing.t.and.disshare.x %>%
    distinct(family.x, family.y, annotation, category) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(family.x, family.y)), collapse = "&", sep = "")) %>% 
    group_by(annotation, category)  %>%
    summarise(num.mosaic.pairs = n_distinct(pair)) %>%
    ungroup() %>%
    left_join(all.family.pairs.sharing.t) %>%
    mutate(prob.mosaicism.if.sharing.fold.and.multi.x.family.measure = num.mosaic.pairs/num.pairs) %>%
    distinct(annotation, category, prob.mosaicism.if.sharing.fold.and.multi.x.family.measure)
  return(mosaicism.per.family)
}

Plot.Domain.Combinations = function(tile.data.object, text.size.axis = 9, text.size.panel = 12, interactive.plot = FALSE) {
  domain.position.plots = list()
  domain.legend.plots = list()
  tile.data.object = tile.data.object %>%
    # keep only the ones with multidomain
    group_by(qname) %>%
    mutate(num.domains = n_distinct(domain[domain != "undetected"])) %>%
    ungroup() %>%
    filter(num.domains > 1)
  
  for (this.annotation in unique(tile.data.object$annotation) %>% sort()) {
    this.tile.data.raw = tile.data.object %>% 
      filter(annotation == this.annotation)  %>% 
      mutate(annotation = this.annotation) %>%
      arrange(row) %>%
      group_by(qname) %>%
      arrange(pos.start, pos.end) %>% 
      mutate(domain.combination = paste(unique(domain), collapse = " | ", sep = ""),
             domain.combination2 = paste(unique(domain[domain != "undetected" & domain != "multiple domains"]), collapse = " | ", sep = "")) %>%
      ungroup() %>%
      arrange(qname) 
    
    
    num.represeq.per.combination =  this.tile.data.raw %>%
      distinct(domain.combination2, qname, n.prot.in.reprseq) %>%
      group_by(domain.combination2) %>%
      summarise(num.prot.in.combination = sum(n.prot.in.reprseq)) %>%
      ungroup()
      
    this.tile.data = this.tile.data.raw %>%
      # select only one qname per domain combination
      group_by(domain.combination2) %>%
      arrange(pos.start) %>%
      mutate(qname1 = qname[1],
             num.reprseq.this.domain.combination = n_distinct(qname)) %>%
      ungroup() %>%
      left_join(num.represeq.per.combination) %>%
      filter(qname == qname1) %>%
      rowwise() %>%
      mutate(domain_new = if_else(
        domain_name == "undetected" | domain_name == "multiple domains",
        domain_name,
        #paste(domain_name, " (X: ", x_name, " [nr  ",x.index ,"]", ")", collapse = "", sep = ""))) %>%
        #paste(x.index, ": ", domain_name, collapse = "", sep = ""))) %>%
        paste(domain, ": ", domain_name,  collapse = "", sep = "")),
        qname = paste(qname, "(", num.prot.in.combination, ")",collapse = "", sep = "")) %>%
      ungroup() #%>%
    #rowwise() %>%
    #mutate(domain_new = gsub(" ", "\n", domain_new))
    
    
    # create a color pallete where domains from common X are different shades of the same color
    this.x = this.tile.data %>% filter(!(domain %in% c("undetected", 
                                                       "multiple domains"))) %>% distinct(x.index) %>% pull(x.index)
    #x.colors = randomcoloR::distinctColorPalette(length(this.x))
    x.colors = colorRampPalette(brewer.pal(12, "Paired"))(length(this.x))
    x.colors = data.frame(x.index = this.x, color = x.colors)
    all.colors = data.frame(domain = character(0), color = character(0))
    for (this.x.index in this.x) {
      all.domains = this.tile.data %>% filter(x.index == this.x.index) %>% distinct(domain_new) %>% 
        pull(domain_new)
      this.color = x.colors %>% filter(x.index == this.x.index) %>% pull(color)
      this.x.colors.pallete =  colorRampPalette(this.color)
      this.colors = this.x.colors.pallete(length(all.domains))
      all.colors = rbind(all.colors, data.frame(domain = all.domains, color = this.colors))
      print(length(all.domains))
    }
    
    all.colors = rbind(all.colors %>% arrange(domain), data.frame(domain = c("undetected", "multiple domains"), 
                                                                  color = c( "#FFFFFF",  "#000000"), stringsAsFactors = FALSE))
    
    colm = all.colors$color
    names(colm) = all.colors$domain
    
    
    this.tile.data.row.numbers = this.tile.data %>%
      arrange(domain.combination2, qname) %>%
      distinct(domain.combination2, qname) %>%
      mutate(dummy = 1) %>%
      group_by(dummy) %>%
      mutate(row.plus = 1:n()) %>%
      mutate(row = row.plus - 0.5)
    
    this.tile.data = this.tile.data %>%
      select(-row, -row.plus) %>%
      left_join(this.tile.data.row.numbers)
    
    
    legend.dummy.plot.data = this.tile.data %>% 
      distinct(x.index, domain_new) %>%
      arrange(x.index)
    dummy.plot = ggplot(legend.dummy.plot.data) + 
      geom_col(aes(x = domain_new, y = 1, fill = domain_new)) +
      scale_fill_manual(values = colm[which(names(colm) %in% unique(legend.dummy.plot.data$domain_new))], 
                        name = "Domain",
                        guide = guide_legend(ncol = 1)) + 
      theme(legend.margin=margin(c(0,0,0,0)),
            legend.key.size = unit(0.45, "cm"))
    leg <- ggpubr::get_legend(dummy.plot, position = "right")
    ggleg = ggpubr::as_ggplot(leg)
    #domain.legend.plots[[this.annotation]] =  ggleg

    num.rows =  this.tile.data %>% distinct(qname) %>% nrow()
    gg_rect = Show.Domain.Position.Within.Data(
      tile.data = this.tile.data %>% select(-domain) %>% rename(domain = domain_new),
      colormap = colm, 
      my.theme = theme.no.verical, 
      color.frames = FALSE,
      interactive.plot = interactive.plot,
      legend.position = "none") +
      #scale_y_discrete(breaks =this.tile.data$qname, labels = this.tile.data$family) +
      
      theme(strip.text.x = element_text(size = text.size.panel), 
            #legend.box="vertical", legend.margin=margin(),
            axis.title.x = element_text(size = text.size.axis),
            axis.title.y = element_text(size = text.size.axis),
            #axis.text.x = element_blank(),
            #axis.text.y = element_blank(),
            axis.text.y = element_text(size = text.size.axis, vjust = 1.1 - (num.rows - 20)/20),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_rect(colour = "gray", fill=NA, size=1)
      ) 
    #if (interactive.plot) {gg_rect = girafe(ggobj=gg_rect)}
    domain.position.plots[[this.annotation]] = plot_grid(gg_rect, ggleg,  ncol = 2,rel_widths = c(0.65, 0.35))
    
  }
  return(list(domain.position.plots = domain.position.plots, gg_rect = gg_rect))
}



Get.Domain.coverage.In.Recent.HGT = function(recent_HGT_pairs_raw, hhr.table, surely.annotated.proteins.from.included.categories, ecod.domains.hits) {
  qnames.in.recent.HGT = unique(c(recent_HGT_pairs_raw$qname, recent_HGT_pairs_raw$sname))
  hhr.table.filtered = hhr.table %>% 
    filter(qname %in% qnames.in.recent.HGT | sname %in% qnames.in.recent.HGT) %>%
    select(qname, sname, qstart, qend, sstart, send, prob)
  
  domains.in.qnames.in.recent.HGT = ecod.domains.hits %>% 
    filter(qname %in% qnames.in.recent.HGT) %>%
    select(qname, domain.hit.qstart = qstart, domain.hit.qend =qend, domain.hit.domain.start = sstart, domain.hit.domain.end = send, dmain.length = slength, t_id, prob) %>%
    # we will choose one t_name per qname, this is a simplification as we may potentially have twice the same t name in one qname
    group_by(qname, t_id) %>%
    arrange(desc(prob)) %>%
    mutate(order = row_number()) %>%
    ungroup() %>%
    filter(order == 1)
  
  
  recent_HGT_pairs_raw.data.1 = recent_HGT_pairs_raw %>%
    select(qname, sname) %>%
    left_join(hhr.table.filtered, by= c('qname', 'sname'))
  
  recent_HGT_pairs_raw.data.2 = recent_HGT_pairs_raw %>%
    select(sname = qname, qname = sname) %>%
    left_join(hhr.table.filtered, by= c('qname', 'sname'))
  
  # make sure every protein in a pair is at least one on a query position
  recent_HGT_pairs_raw.data = rbind(recent_HGT_pairs_raw.data.1, recent_HGT_pairs_raw.data.2) 
  nrow.missing = nrow(recent_HGT_pairs_raw.data %>% filter(is.na(prob))) 
  if ( nrow.missing > 1) {
    print(paste0(nrow.missing, " some pairs are missing from hhr table"))
    recent_HGT_pairs_raw.data = recent_HGT_pairs_raw.data %>% filter(!is.na(prob))
  }
  
  recent_HGT_pairs_raw.domain.data.raw = recent_HGT_pairs_raw.data %>%
    distinct(qname, sname, qstart, qend) %>%
    # we may have multiple domains in one qname
    left_join(domains.in.qnames.in.recent.HGT, by = "qname") %>%
    rowwise() %>%
    mutate(fragment.length = qend - qstart + 1,
           domain.in.frament.length = min(domain.hit.qend, qend) - max(qstart, domain.hit.qstart) + 1) %>%
    mutate(detected.domain.within.shared.fragment = domain.in.frament.length > 0,
           domain.in.frament.l = max(domain.in.frament.length, 0),
           prop.domain.in.fragment = domain.in.frament.l/dmain.length,
           prop.fragment.covered.by.domain = domain.in.frament.length/fragment.length) %>%
    ungroup()

  recent_HGT_pairs_raw.domain.data = recent_HGT_pairs_raw.domain.data.raw %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(qname, sname)), collapse = " AND ", sep = "")) %>% 
    filter(detected.domain.within.shared.fragment) %>%
    group_by(pair, t_id) %>%
    summarise(
      # PREVIOUSLY WE USED MEAN!!
      prop.domain.in.fragment = max(prop.domain.in.fragment, na.rm = TRUE),
      prop.fragment.covered.by.domain = max(prop.fragment.covered.by.domain, na.rm = TRUE),
      n = n()) %>%
    ungroup()
    
    # make it symmetric
    recent_HGT_pairs_raw.domain.data1 = recent_HGT_pairs_raw.domain.data %>% tidyr::separate(col = "pair", into = c("qname", "sname"), sep = " AND ")
    recent_HGT_pairs_raw.domain.data2 = recent_HGT_pairs_raw.domain.data %>% tidyr::separate(col = "pair", into = c("sname", "qname"), sep = " AND ")
  
    recent_HGT_pairs_raw.domain.data.to.return = rbind(recent_HGT_pairs_raw.domain.data1, recent_HGT_pairs_raw.domain.data2) %>%
    left_join(surely.annotated.proteins.from.included.categories %>% select(qname, qannot = annotation))  %>%
    left_join(surely.annotated.proteins.from.included.categories %>% select(sname = qname, sannot = annotation))
  return(recent_HGT_pairs_raw.domain.data.to.return)
}
  
  
  
Get.Recent.HGT.Domains = function(recent_HGT_pairs_raw, hhr.table, surely.annotated.proteins.from.included.categories, ecod.domains.hits) {
  recent_HGT_pairs_raw.domain.data = Get.Domain.coverage.In.Recent.HGT(recent_HGT_pairs_raw, hhr.table, surely.annotated.proteins.from.included.categories, ecod.domains.hits)
  
  # there should be max 2 hits per pair
  if (nrow( recent_HGT_pairs_raw.domain.data %>% filter(n > 2)) > 0) {print("duplicated entries")}
  
  recent_HGT_pairs_raw.domain.data.to.plot.all = recent_HGT_pairs_raw.domain.data %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(qname, sname)), collapse = "&", sep = "")) %>% 
    group_by(pair, t_id) %>%
    # note that they may actually differ a lot if we take difernt f_id fir both sides: they may have different length
    summarise(av.prop.domain.in.fragment = mean(prop.domain.in.fragment),
              av.prop.fragment.covered.by.domain = mean(prop.fragment.covered.by.domain)) %>%
    group_by(t_id) %>%
    summarise(n_pairs = n_distinct(pair), 
              av.prop.domain.in.fragment = round(mean(av.prop.domain.in.fragment),2),
              av.prop.fragment.covered.by.domain = round(mean(av.prop.fragment.covered.by.domain),2)) %>%
    arrange(desc(n_pairs)) 
  
  
  
  recent_HGT_pairs_raw.domain.data.to.plot.all$h.index <- sapply(1:nrow(recent_HGT_pairs_raw.domain.data.to.plot.all), function(index){
    this.topology <- recent_HGT_pairs_raw.domain.data.to.plot.all$t_id[index]
    this.topology.h.index <- strsplit(this.topology, ".", fixed = T)[[1]]
    this.topology.h.index <- this.topology.h.index[1:2]
    this.topology.h.index <- as.character(paste(this.topology.h.index, collapse = "."))
  })
  
  recent_HGT_pairs_raw.domain.data.to.plot.all = recent_HGT_pairs_raw.domain.data.to.plot.all %>%
    left_join(ecod_id_to_names_map %>% filter(level == "H") %>% select(h.index = id, h_name = name)) %>%
    left_join(ecod_id_to_names_map %>% filter(level == "T") %>% select(t_id = id, t_name = name)) 
  
  return(recent_HGT_pairs_raw.domain.data.to.plot.all)
}




Plot.Domains.Within.Mosaic.Proteins = function(domains.within.proteins.domain.mosaic, colors = NULL) {
  
  if (is.null(colors)) {
    colors = domains.within.proteins.domain.mosaic %>% distinct(t_id) %>% mutate(color = "gray")
  }
  
  domains.within.proteins.with.significant.higher.domain.mosaic = domains.within.proteins.domain.mosaic %>%
    left_join(colors, by = "t_id") %>%
    filter(corrected.pval < 0.05) %>%
    select(t_id, t_name, h_name, mosaic,odds.ratio, corrected.pval, fillcolor = color) %>%
    arrange(desc(odds.ratio))
  #write.xlsx(x = domains.within.proteins.with.significant.higher.domain.mosaic, 
  #           file = sprintf("%stables/Supplementary_table_mosaic_vs_ninmosaic-domain.xlsx", OUTPUT.FIGURES.PATH), 
  #           sheetName = "Domains in mosaic vs non-mosaic proteins", 
  #  col.names = TRUE, row.names = TRUE, append = TRUE)
  
  
  included.ecods = domains.within.proteins.with.significant.higher.domain.mosaic %>% 
    filter(corrected.pval < 0.05) %>%
    arrange(desc(odds.ratio)) %>%
    head(20) %>%
    #filter(freq.q >= 0.025) %>%
    pull(t_id) %>% unique()
  
  ecod.table <- domains.within.proteins.with.significant.higher.domain.mosaic %>% 
    filter(t_id %in% included.ecods )%>%
    arrange(mosaic, odds.ratio)
  
  ecod.table$t_name = factor(ecod.table$t_name, 
                             levels = ecod.table %>% filter(mosaic == "Mosaic") %>% arrange(odds.ratio) %>% pull(t_name) %>% unique())
  #ecod.table$mosaic = factor(ecod.table$mosaic, levels = c("Not mosaic", "Mosaic"))
  
  display.h.name <- ecod.table$h_name
  display.h.name[ecod.table$t_name == ecod.table$h_name] <- "same"
  display.t.name <- sprintf("%s\n(H: %s)", ecod.table$t_name, display.h.name)
  
  
  p.odds = ggplot(ecod.table %>% distinct(t_name, odds.ratio, fillcolor), aes(x=t_name, y = odds.ratio, fill = fillcolor)) +
    geom_bar(stat="identity", width = 0.8, col = "black") + 
    scale_x_discrete(breaks = ecod.table$t_name, labels = display.t.name) +
    scale_fill_identity() +
    theme_minimal() + 
    theme(legend.text=element_text(size=4)) + 
    labs(x = "", y = "Odds ratio", fill = "category") + 
    coord_flip() +
    guides(fill = "none")
  return(p.odds)
}


Get_Domain_Odss_In_Fold_Types = function(all.fold.types, ecod_id_to_names_map) {
  all.fold.types.stats = all.fold.types %>%
    group_by(fold.type.mosaic) %>%
    mutate(num.this.mosaic.fold.type = n_distinct(fold.type)) %>%
    group_by(t_id) %>%
    mutate(num.fold.types.with.this.t = n_distinct(fold.type)) %>%
    filter(num.fold.types.with.this.t >= 1) %>%
    group_by(t_id, fold.type.mosaic, num.this.mosaic.fold.type) %>%
    summarise(num.fold.types.with.this.t = n_distinct(fold.type)) %>%
    ungroup() %>%
    tidyr::complete(t_id, fold.type.mosaic, 
                    fill = list(num.fold.types.with.this.t = 0, num.this.mosaic.fold.type = 0)) %>%
    mutate(num.fold.types.with.other.t  = num.this.mosaic.fold.type - num.fold.types.with.this.t) %>% 
    ungroup()
  
  all.t.in.fold.types = all.fold.types.stats %>% distinct(t_id)
  all.t.in.fold.types$pval = NA
  all.t.in.fold.types$odds.ratio = NA
  for (this.t.id in unique(all.fold.types.stats$t_id)) {
    this.t.stats = all.fold.types.stats %>% filter(t_id == this.t.id)
    M1 = this.t.stats %>% distinct(fold.type.mosaic, num.fold.types.with.this.t, num.fold.types.with.other.t) %>% as.data.table()
    #domains.within.proteins$pval[which(domains.within.proteins$t_id == this.t.id)] = chisq.test(M1)$p.value
    all.t.in.fold.types$pval[which(all.t.in.fold.types$t_id == this.t.id)] = fisher.test(M1, alternative = "less")$p.value
    all.t.in.fold.types$odds.ratio[which(all.t.in.fold.types$t_id == this.t.id)] = (M1$num.fold.types.with.this.t[2]/M1$num.fold.types.with.other.t[2])/(M1$num.fold.types.with.this.t[1]/M1$num.fold.types.with.other.t[1])
  }
  
  all.t.in.fold.types$h.index <- sapply(1:nrow(all.t.in.fold.types), function(index){
    this.topology <- all.t.in.fold.types$t_id[index]
    this.topology.h.index <- strsplit(this.topology, ".", fixed = T)[[1]]
    this.topology.h.index <- this.topology.h.index[-length(this.topology.h.index)]
    this.topology.h.index <- as.character(paste(this.topology.h.index, collapse = "."))
  })
  
  all.t.in.fold.types = all.t.in.fold.types %>% 
    left_join(ecod_id_to_names_map %>% filter(level == "H") %>% select(h.index = id, h_name = name)) %>%
    left_join(ecod_id_to_names_map %>% filter(level == "T") %>% select(t_id = id, t_name = name)) %>%
    mutate(corrected.pval = pval/nrow(all.t.in.fold.types)) %>%
    mutate(mosaic = "Mosaic") %>%
    #filter(pval < 0.05/nrow(all.t.in.fold.types)) %>%
    arrange(desc(odds.ratio))
  return(all.t.in.fold.types)
}



Get.Mosaicism.Odds.Ratio.Data = function(all.proteins, annotations, diversity.data) {
  all.proteins.with.any.relaxed.function = all.proteins %>%
    distinct(qname, mosaic) %>%
    left_join(annotations) %>%
    filter(!is.na(annotation.index)) %>%
    group_by(annotation.index) %>%
    mutate(num.proteins.in.this.annotation = n_distinct(qname)) %>%
    ungroup() 
  
  num.proteins.per.osaicism.data = all.proteins.with.any.relaxed.function %>%
    group_by(mosaic) %>%
    summarise(num.proteins = n_distinct(qname)) %>%
    ungroup()
  
  annot.stats.per.category.and.mosaicism = all.proteins.with.any.relaxed.function %>%
    inner_join(diversity.data) %>%
    group_by(annotation.index, mosaic, annotation, category) %>%
    summarise(num.proteins.with.this.annot = n_distinct(qname)) %>%
    ungroup() %>%
    left_join(num.proteins.per.osaicism.data) %>%
    mutate(num.proteins.with.different.annot = num.proteins - num.proteins.with.this.annot,
           prop.proteins.has.this.annot = num.proteins.with.this.annot/num.proteins) %>%
    tidyr::complete(annotation.index, mosaic, fill = list(num.proteins = 0, num.proteins.with.this.annot = 0, num.proteins.with.different.annot = 0, prop.proteins.has.this.annot = NA))
  
  
  
  annot.stats.per.category.and.mosaicism$pval = NA
  annot.stats.per.category.and.mosaicism$odds.ratio = NA
  for (this.annotation.index in unique(annot.stats.per.category.and.mosaicism$annotation.index)) {
    this.annot.stats = annot.stats.per.category.and.mosaicism %>% filter(annotation.index == this.annotation.index)
    M1 = this.annot.stats %>% select(mosaic, num.proteins.with.this.annot, num.proteins.with.different.annot) %>% select(num.proteins.with.this.annot, num.proteins.with.different.annot) %>% as.data.table()
    #  annot.stats.per.category.and.mosaicism$pval[which(annot.stats.per.category.and.mosaicism$annotation.index == this.annotation.index)] = chisq.test(M1)$p.value
    annot.stats.per.category.and.mosaicism$pval[which(annot.stats.per.category.and.mosaicism$annotation.index == this.annotation.index)] = fisher.test(M1)$p.value
    annot.stats.per.category.and.mosaicism$odds.ratio[which(annot.stats.per.category.and.mosaicism$annotation.index == this.annotation.index)] = (M1$num.proteins.with.this.annot[2]/M1$num.proteins.with.different.annot[2])/(M1$num.proteins.with.this.annot[1]/M1$num.proteins.with.different.annot[1])
  }
  num.comparisons = n_distinct(annot.stats.per.category.and.mosaicism$annotation.index)
  annot.stats.per.category.and.mosaicism = annot.stats.per.category.and.mosaicism %>%
    mutate(corrected.pval = pval*num.comparisons)
  return(annot.stats.per.category.and.mosaicism)
}


GetDomainsWithinMosaicProteins = function(all.proteins, families, ecod.domains.hits, ecod_id_to_names_map) {
  
  domains.within.proteins = all.proteins %>%
    # filter(has.hit.to.anything.in.the.data) %>%
    left_join(families) %>%
    inner_join(ecod.domains.hits, by = "qname") %>%
    group_by(mosaic) %>%
    mutate(total.num.families.with.any.ecod.hit = n_distinct(family),
           total.num.proteins.with.any.ecod.hit = n_distinct(qname)) %>%
    ungroup() %>%
    distinct(qname, family, t_id, t_name, mosaic, total.num.families.with.any.ecod.hit, total.num.proteins.with.any.ecod.hit) %>%
    group_by(t_id, t_name, mosaic, total.num.families.with.any.ecod.hit, total.num.proteins.with.any.ecod.hit) %>%
    summarise(num.families = n_distinct(family),
              num.proteins.with.this.ecod = n_distinct(qname)) %>%
    ungroup() %>%
    mutate(freq = num.families/total.num.families.with.any.ecod.hit,
           freq.q = num.proteins.with.this.ecod/total.num.proteins.with.any.ecod.hit,
           num.proteins.other.ecod.hits = total.num.proteins.with.any.ecod.hit - num.proteins.with.this.ecod) %>%
    tidyr::complete(t_id, mosaic, fill = list(num.proteins.with.this.ecod = 0, num.proteins.other.ecod.hits = 0))
  
  
  domains.within.proteins$h.index <- sapply(1:nrow(domains.within.proteins), function(index){
    this.topology <- domains.within.proteins$t_id[index]
    this.topology.h.index <- strsplit(this.topology, ".", fixed = T)[[1]]
    this.topology.h.index <- this.topology.h.index[-length(this.topology.h.index)]
    this.topology.h.index <- as.character(paste(this.topology.h.index, collapse = "."))
  })
  
  
  
  domains.within.proteins$pval = NA
  domains.within.proteins$odds.ratio = NA
  for (this.t.id in unique(domains.within.proteins$t_id)) {
    this.domain.stats = domains.within.proteins %>% filter(t_id == this.t.id)
    M1 = this.domain.stats %>% select(mosaic, num.proteins.with.this.ecod, num.proteins.other.ecod.hits) %>% select(num.proteins.with.this.ecod, num.proteins.other.ecod.hits) %>% as.data.table()
    #domains.within.proteins$pval[which(domains.within.proteins$t_id == this.t.id)] = chisq.test(M1)$p.value
    domains.within.proteins$pval[which(domains.within.proteins$t_id == this.t.id)] = fisher.test(M1)$p.value
    domains.within.proteins$odds.ratio[which(domains.within.proteins$t_id == this.t.id)] = (M1$num.proteins.with.this.ecod[2]/M1$num.proteins.other.ecod.hits[2])/(M1$num.proteins.with.this.ecod[1]/M1$num.proteins.other.ecod.hits[1])
  }
  
  num.comparisons = n_distinct(domains.within.proteins$t_id)
  domains.within.proteins = domains.within.proteins %>%
    mutate(corrected.pval = pval*num.comparisons)
  
  domains.within.proteins = domains.within.proteins %>%
    left_join(ecod_id_to_names_map %>% filter(level == "H") %>% select(h.index = id, h_name = name)) %>%
    mutate(mosaic = if_else(mosaic, "Mosaic", "Not mosaic"))
  return(domains.within.proteins)
}

Plot.Mosaicism.Table = function(diversity.top, PHROG.COLOR.MAP) {
  p <- ggplot(diversity.top, aes(y = annotation, x = value, fill = category)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    scale_fill_manual(values = PHROG.COLOR.MAP) + theme_minimal() + 
    theme(strip.text.y = element_text(face = "bold", size = 7, angle = 360)) + 
    labs(x = "", y = "") +
    scale_x_continuous(labels=function(x) sprintf("%.2f", x))
  
  p.most.diverse <- p + guides(fill = "none")
  return(p.most.diverse)
}



unique.unclassified.rm = function(x, unknown_words = c("unclassified", "unknown", "unspecified", ""), unknown_word_to_use = "unknown") {
    classified.x.index = which(!(tolower(x) %in% unknown_words))
    if (length(classified.x.index) == 0) {
      return(unknown_word_to_use)
    } else {
      y = x[classified.x.index]
      if (length(unique(y)) == 1) {
        return(unique(y))
      } else {
        return("multi")
      }
    }
}
    
