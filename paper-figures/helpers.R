
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
                               edge.colors,
                               node.color.map,
                               node.size.map = NULL,
                               line.color.variable,
                               vertex.size = 1) {
  if (is.null(node.size.map)) {
    node.size.map = data.frame(node = unique(c(data.for.network$from, data.for.network$to))) %>% mutate(size =vertex.size)
  }
  label.color.map = node.color.map
  data.for.network$edge.var =  data.for.network %>% pull(line.color.variable)
  net_hh = igraph::graph_from_data_frame(data.for.network)# %>% filter(min.cov == 0.3))
  df_edges <- igraph::as_data_frame(net_hh, what = "edges") %>% 
    select(from, to, edge.var)
  
  net_hh = igraph::simplify(net_hh,
                              remove.multiple = TRUE,
                              remove.loops = TRUE)
    
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
    igraph::E(net_annotated)$lty = "solid"
  } else {
    for (i in 1:nrow(edge.colors)) {
      igraph::E(net_annotated)$edge.color[igraph::E(net_annotated)$edge.var == edge.colors$value[i]] <- edge.colors$color[i]
      igraph::E(net_annotated)$lty[igraph::E(net_annotated)$edge.var == edge.colors$value[i]] <- edge.colors$lty[i]
    }
  }
  return(net_annotated)
}


#' Visualie network using plotly package.
#' We need to first convert our network to a data.frame using ggnetwork package
#' @net_annotated network with annotation attribute
#' @colors a mapping of annotations and colors
VisualiseGGnetwork = function(net_annotated, vertex.size = 0.5, cex = 0.1, label_size = 1, layout = igraph::nicely(), main = "") {
  data = ggnetwork(x = igraph::simplify(net_annotated, remove.multiple = T, remove.loops = T), 
                   layout = layout,
                   scale = TRUE, stringsAsFactors = FALSE, arrow.gap = 0)
  
  
  if (!("size" %in% names(data))) {
    data$size = vertex.size
  }
  if (!("color" %in% names(data))) {
    data$color = "gray"
  }
  
  g=ggplot(data) +
    geom_edges(aes(x = x, y = y, xend=xend, yend=yend), color = "grey50", size = cex) +
    geom_nodes(aes(x = x, y = y, color = color, label1 = name, size = size)) +
    geom_text(aes(x = x, y = y, color = color, label = name), size = label_size, nudge_y = 0.025) +
    guides(size = "none") +
    scale_color_identity() +
    xlim(c(-0.05, 1.05)) +
    theme_blank() +
    ggtitle(main)
  
  return(g)
}


Get.Ecod.Annotation.Mapping.Numeric = function(distinct_domains, ecod.domain.descriptions) {

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
  ecod.annotation.map$X.uid =  sapply(1:nrow(ecod.annotation.map),
                                      function(row.index){Get.uid.From.ECOD.subject.id(ecod.annotation.map$domain[row.index])})
  ecod.annotation.map$ecod_domain_id = sapply(1:nrow(ecod.annotation.map),
                                              function(row.index){Get.ecod.domain.id.From.ECOD.subject.id(ecod.annotation.map$domain[row.index])})
  ecod.annotation.map = ecod.annotation.map %>%
    left_join(ecod.domain.descriptions.full.names, by = c("X.uid", "ecod_domain_id"))
  return(ecod.annotation.map)
}

Get.Domain.Positions.Data = function(annotated.ecod.domains) {
  annotated.ecod.domains = annotated.ecod.domains %>% 
    select(qname, annotation, qlength, qstart, qend, domain, family) %>%
    arrange(family)
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


Show.Domain.Position.Within.Data = function(tile.data, colormap, my.theme = theme_bw()) {
  if (nrow(tile.data) == 0) {
    gg_rect = ggplot() + geom_blank()
  } else {
    tile.data = tile.data %>%
      mutate(tooltip = paste0("Domain annotation:", domain_name, "\n", "seq: ", qname, "\n", "domain.id: ", domain)) 
    
    gg_rect = ggplot(tile.data) +
      geom_blank(aes(x=pos.start, y = qname)) +
      scale_x_continuous(name="position") +
      ylab("") +
      geom_rect_interactive(mapping = aes(xmin = pos.start, xmax = pos.end,
                                          ymin = row, ymax = row.plus, fill = domain, color=domain,
                                          tooltip = tooltip), alpha=1) +
      # leave in the egend only the colors that we can actually see
      scale_fill_manual(values = colormap[which(names(colormap) %in% unique(tile.data$domain))], name = "Domain") + 
      scale_color_manual(values = colormap[which(names(colormap) %in% unique(tile.data$domain))], name = "Domain")  +
      facet_grid(. ~annotation) + 
      my.theme +
      guides(
        color = guide_legend(show = FALSE) 
      )
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

Calculate.Annotation.Tradeoff = function(hhr.phrogs.annotated, max.eval.range, min.cov.range, min.prob.range) {
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


