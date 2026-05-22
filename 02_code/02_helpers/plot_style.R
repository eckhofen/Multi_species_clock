library(ggplot2)
library(patchwork)

# plot settings
theme_custom <- theme_minimal() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            # panel.background = element_rect(color = "black", linewidth = 1),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12), 
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(0.1, "cm"), 
            axis.title = element_text(size = 10), 
            legend.ticks = element_blank())

theme_set(theme_custom)

save_plot <- function(filename, plot = last_plot(), width = 10, height = 8, dpi = 300, format = c("png", "pdf"), ...) {
  formats <- c("png", "pdf", "jpg")
  selected_formats <- formats[formats %in% format]
  filename <- tools::file_path_sans_ext(filename)
  for (fmt in selected_formats) {
    ggsave(paste0(filename, ".", fmt), plot, width = width, height = height, dpi = dpi, ...)
  }
}

