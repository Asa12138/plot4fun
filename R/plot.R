#' Word cloud plot
#'
#' @param str_vector string vector
#' @param ignore_words ignore_words
#' @param topN topN, 50
#'
#' @export
#' @return a htmlwidget
#' @examples
#' \donttest{
#' data(otutab, package = "pcutils")
#' if (requireNamespace("wordcloud2")) {
#'   my_wordcloud(taxonomy$Genus)
#' }
#' }
my_wordcloud <- function(str_vector,
                         ignore_words = "Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome|Unassig", topN = 50) {
  lib_ps("wordcloud2", library = FALSE)
  str_vector <- str_vector[!grepl(ignore_words, str_vector)]
  sort(table(str_vector), decreasing = TRUE)[1:topN] %>%
    as.data.frame() %>%
    stats::na.omit() %>%
    wordcloud2::wordcloud2(., size = .7)
}

#' Give you a rose
#'
#' @param color "skyblue3"
#'
#' @return plot
#' @export
#' @references <https://mp.weixin.qq.com/s/W-BYPR3UXL120XWpTmN3rA>
give_you_a_rose <- function(color = "red3") {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mar = rep(1, 4))

  lib_ps("plot3D", library = FALSE)
  # 生成绘图数据
  x <- seq(0, 24) / 24
  t <- seq(0, 575, by = 0.5) / 575 * 20 * pi + 4 * pi
  grid <- expand.grid(x = x, t = t)
  x <- matrix(grid$x, ncol = 25, byrow = TRUE)
  t <- matrix(grid$t, ncol = 25, byrow = TRUE)
  p <- (pi / 2) * exp(-t / (8 * pi))
  change <- sin(15 * t) / 150
  u <- 1 - (1 - (3.6 * t) %% (2 * pi) / pi)^4 / 2 + change
  y <- 2 * (x^2 - x)^2 * sin(p)
  r <- u * (x * sin(p) + y * cos(p))
  # 绘图
  plot3D::persp3D(
    x = r * cos(t), y = r * sin(t), z = u * (x * cos(p) - y * sin(p)),
    main = "To you",
    # xlim=c(-0.5,0.5),ylim=c(-0.5,0.5),zlim=c(0,1),
    xlab = "Love youself",
    ylab = "Love youself",
    zlab = "Love youself",
    col = grDevices::colorRampPalette(c("#e4e9f6", color))(100),
    border = "grey85",
    lwd = 0.1,
    facets = TRUE,
    colkey = FALSE,
    bty = "b2",
    theta = -60, phi = 45
  )
  message("give you a rose \ud83c\udf39.")
}


#' Plot a DNA double helix
#'
#' @param col_DNA col_DNA, "#377EB8"
#' @param col_ATCG col_ATCG, c("#7FC97F","#FB8072","#FFFFB3","#A6CEE3")
#' @param DNA_length DNA_length, 2
#' @return ggplot
#' @export
#' @references <https://github.com/SherryDong/create_plot_by_R_base>
#' @examples
#' DNA_plot()
DNA_plot <- function(col_DNA = "#377EB8", col_ATCG = c("#7FC97F", "#FB8072", "#FFFFB3", "#A6CEE3"),
                     DNA_length = 2) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphics::par(pin = c(1.5 * DNA_length, 1.5))

  DNA_length <- 2 * DNA_length ## the code only applies when DNA_length%%2==0, if DNA_length%%2==1, need to modify

  x <- seq(-DNA_length * pi / 2, DNA_length * pi / 2, length.out = 1000) ##
  y1 <- cos(x) ## backbone up
  y2 <- cos(x + pi) ## backbone down
  # get the position of nucleotides
  xx <- seq(DNA_length * pi / 2, -DNA_length * pi / 2, length.out = DNA_length * 5 + 1)
  xx <- xx + (xx[2] - xx[1]) / 2
  # remove the first and the lines in the boundary region
  xx <- setdiff(xx, c(xx[c(1:DNA_length) * 5 - 2], min(xx)))
  plot(y1 ~ x, pch = 16, type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "", bty = "n", col = "white")
  for (i in 1:length(xx)) {
    ybottom <- cos(xx[i]) # ybottom position
    ytop <- cos(xx[i] + pi) # yup position
    rr <- sample(1:4, 1) ## ATCG, random select one pair
    if (rr == 1) {
      graphics::segments(y0 = ybottom, y1 = 0, x0 = xx[i], x1 = xx[i], col = col_ATCG[1], lwd = 4)
      graphics::segments(y0 = 0, y1 = ytop, x0 = xx[i], x1 = xx[i], col = col_ATCG[2], lwd = 4)
    }
    if (rr == 2) {
      graphics::segments(y0 = ybottom, y1 = 0, x0 = xx[i], x1 = xx[i], col = col_ATCG[2], lwd = 4)
      graphics::segments(y0 = 0, y1 = ytop, x0 = xx[i], x1 = xx[i], col = col_ATCG[1], lwd = 4)
    }
    if (rr == 3) {
      graphics::segments(y0 = ybottom, y1 = 0, x0 = xx[i], x1 = xx[i], col = col_ATCG[3], lwd = 4) ## C-G
      graphics::segments(y0 = 0, y1 = ytop, x0 = xx[i], x1 = xx[i], col = col_ATCG[4], lwd = 4)
    }
    if (rr == 4) {
      graphics::segments(y0 = ybottom, y1 = 0, x0 = xx[i], x1 = xx[i], col = col_ATCG[4], lwd = 4) ## G-C
      graphics::segments(y0 = 0, y1 = ytop, x0 = xx[i], x1 = xx[i], col = col_ATCG[3], lwd = 4)
    }
  }
  graphics::lines(y1 ~ x, pch = 16, lwd = 8, col = col_DNA)
  graphics::lines(y2 ~ x, pch = 16, lwd = 8, col = col_DNA)
}

#' Draw a Chunlian (Spring Festival couplet) using ggplot2
#'
#' @param words A character vector containing three strings for the three lines of the couplet
#' @param bg_size Size of the points in geom_point, 20
#' @param bg_shape Shape of the points in geom_point (21~25), 22 or 23 are very good.
#' @param bg_fill Fill color of the points in geom_point
#' @param font_file font file, e.g XX.ttf, XX.ttc
#' @param download_dir download_dir for font_file
#' @param text_size Size of the text in geom_text, 10
#' @param text_params parameters parse to geom_text
#'
#' @return A ggplot object representing the Chunlian
#' @export
chunlian <- function(words = NULL, bg_size = 20, bg_shape = 22, bg_fill = "red2", text_size = 10, text_params = list(), font_file = NULL, download_dir = "plot4fun_temp") {
  x <- y <- label <- NULL
  if (identical(words, 1)) {
    words <- c("\u6295\u5565\u4e2d\u5565", "SCI\u5929\u5929\u6709\u4e00\u533a", "CNS\u6708\u6708\u6709\u5c01\u9762")
  } else if (identical(words, 2)) words <- c("\u79d1\u7814\u987a\u5229", "\u6570\u636e\u5206\u6790\u597d\u5230\u7206", "\u6587\u7ae0\u6295\u54ea\u54ea\u90fd\u8981")

  lib_ps("sysfonts", "showtext", library = FALSE)

  if (is.null(font_file)) {
    font_file <- file.path(download_dir, "SanJiChunLianZiTiJian.ttf")
  }
  if (!file.exists(font_file)) {
    chunlian_font_url <- "https://asa12138.github.io/FileList/SanJiChunLianZiTiJian.ttf"
    download2(chunlian_font_url, font_file)
  }

  if (!file.exists(font_file)) stop("font_file don't exsit.")

  showtext::showtext_auto()
  # Add the font to showtext
  sysfonts::font_add("chunlian", font_file)

  words <- words[1:3]
  words[is.na(words)] <- ""

  hengpi_df <- shanglian_df <- xialian_df <- data.frame()
  hengpi <- strsplit(words[1], "")[[1]]
  if (length(hengpi) > 0) hengpi_df <- data.frame(y = 1, x = seq_along(hengpi), label = hengpi)

  shanglian <- strsplit(words[2], "")[[1]]
  xialian <- strsplit(words[3], "")[[1]]
  if (length(shanglian) > 0) shanglian_df <- data.frame(x = 0, y = -seq_along(shanglian) + 0.5, label = shanglian)
  if (length(xialian) > 0) xialian_df <- data.frame(x = nrow(hengpi_df) + 1, y = -seq_along(xialian) + 0.5, label = xialian)

  dat <- rbind(hengpi_df, shanglian_df, xialian_df)

  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(size = bg_size, shape = bg_shape, fill = bg_fill, color = "NA") +
    do.call(geom_text, update_param(list(
      mapping = aes(label = label),
      size = text_size, family = "chunlian"
    ), text_params)) +
    xlim(range(c(dat$x - 1, dat$x + 1))) +
    ylim(range(c(dat$y - 1, dat$y + 1))) +
    theme_void() +
    coord_fixed()
  p
}

#' Plot the Olympic rings
#'
#' @return ggplot
#' @export
#'
#' @examples
#' Olympic_rings()
Olympic_rings <- function() {
  radius <- x <- y <- color <- start <- end <- NULL

  lib_ps("ggforce", library = FALSE)

  r <- 1
  pensize <- r / 6
  rings_data <- data.frame(
    x = c(-2 * (r + pensize), -(r + pensize), 0, (r + pensize), 2 * (r + pensize)),
    y = c(r, 0, r, 0, r),
    radius = rep(r, 5),
    color = c("#0081C8", "#FCB131", "#000000", "#00A651", "#EE334E")
  )
  tao_data <- data.frame(
    x = c(-(r + pensize), -(r + pensize), (r + pensize), (r + pensize)),
    start = c(0, 5 / 4 * pi, 0, 5 / 4 * pi),
    end = c(1 / 4 * pi, 7 / 4 * pi, 1 / 4 * pi, 7 / 4 * pi),
    color = c("#FCB131", "#FCB131", "#00A651", "#00A651")
  )
  ggplot() +
    ggforce::geom_circle(
      data = rings_data[c(2, 4), ],
      mapping = aes(r = radius, x0 = x, y0 = y, size = I(5), color = color)
    ) +
    ggforce::geom_circle(
      data = rings_data[c(1, 3, 5), ],
      mapping = aes(r = radius, x0 = x, y0 = y, size = I(5), color = color)
    ) +
    ggforce::geom_arc(data = tao_data, mapping = aes(
      x0 = x, y0 = 0, r = r, size = I(5),
      start = start, end = end, color = color
    )) +
    scale_color_identity() +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
}


#' convert a imgage to 01 matrix
#'
#' @param image_file image_file
#' @param size 32
#' @param breaks breaks, default 2
#'
#' @return chr_mat
#' @export
convert_img_to_matrix <- function(image_file, size = 32, breaks = 2) {
  lib_ps("magick", library = FALSE)
  # 读取并返回图像
  image <- magick::image_read(image_file)
  # 将图像转换为灰度并调整大小
  image <- magick::image_convert(image, format = "png", depth = 8)
  image <- magick::image_scale(image, paste0(size, "x", size))

  # 提取像素值并转换为01矩阵
  image_mat <- as.integer(magick::image_data(image)[1, , ])
  image_mat_f <- cut(image_mat, breaks = breaks)
  image_mat <- nlevels(image_mat_f) - as.integer(image_mat_f)
  chr_mat <- matrix(image_mat,
    nrow = size, byrow = TRUE
  )
  class(chr_mat) <- c("chr_mat", class(chr_mat))
  chr_mat
}

#' convert a character to 01 matrix
#'
#' @param char a character
#' @param size 32
#' @param font_file font_file
#' @param picture_dir where to save the temporary picture
#'
#' @return chr_mat
#' @export
#'
#' @examples
#' convert_chr_to_matrix("A")
convert_chr_to_matrix <- function(char, size = 32, font_file = NULL,
                                  picture_dir = tempdir()) {
  # 将汉字渲染成图像的函数
  stopifnot(nchar(char) == 1)
  lib_ps("showtext", "sysfonts", library = FALSE)
  showtext::showtext_auto()
  if (!is.null(font_file)) {
    sysfonts::font_add("new_font", font_file)
  }

  # 使用png图形设备创建图像
  width <- round(256 * size / 16)
  grDevices::png(file = file.path(picture_dir, "temp.png"), width = width, height = width, bg = "white")
  # plot.new()
  # oldpar <- graphics::par(no.readonly = TRUE)
  # on.exit(graphics::par(oldpar))
  # graphics::par(mar = rep(0,4))
  # if(!is.null(font_file))text(0.5, 0.5, char, cex=size, family="new_font")
  # else text(0.5, 0.5, char, cex=size)
  print(ggplot() +
    annotate("text", 0, 0, label = char, size = size * 5) +
    theme_void())
  grDevices::dev.off()

  chr_mat <- convert_img_to_matrix(file.path(picture_dir, "temp.png"), size = size)
  attributes(chr_mat)$name <- char
  return(chr_mat)
}

#' Plot a chr_mat
#'
#' @param colors c("grey","red2")
#' @param x chr_mat object
#' @param ... add
#' @param random add random
#'
#' @return plot
#' @exportS3Method
#' @method plot chr_mat
plot.chr_mat <- function(x, colors = c("grey", "red2"), random = FALSE, ...) {
  d <- x
  if (random) d <- d + stats::rnorm(length(d), 0.2, 0.1)
  lib_ps("reshape2", library = FALSE)
  col <- row <- value <- NULL
  d %>% as.data.frame() -> d
  rownames(d) <- as.character(rownames(d))
  colnames(d) <- as.character(colnames(d))

  rownames(d) -> d$row

  dd <- reshape2::melt(d, id.vars = "row", variable.name = "col")
  dd$row <- factor(dd$row, levels = rev(rownames(d)))
  dd$col <- factor(dd$col, levels = colnames(d))

  p <- ggplot(dd, aes(x = col, y = row, fill = value)) +
    do.call(geom_tile, update_param(list(color = "white"), list(...))) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    xlab(NULL) +
    ylab(NULL) +
    scale_fill_gradientn(colours = colors) +
    coord_fixed()
  p
}

#' make a LED screen
#'
#' @param chars chars
#' @param colors c("grey","red2")
#' @param save_file save_file
#' @param speed pixel speed, default 32
#' @param ... add
#' @param LED_width LED_width
#' @param fps frame per second, 10
#' @param LED_height LED_height, 64
#' @param image_scale image scale, 10
#'
#' @return gif file
#' @export
#'
#' @examples
#' \donttest{
#' if (interactive()) make_LED("SOS!")
#' }
make_LED <- function(chars = "SOS!", save_file = NULL, LED_width = NULL,
                     speed = 32, fps = 10, colors = c("grey", "red2"),
                     LED_height = 32, image_scale = 10, ...) {
  lib_ps("gifski", library = FALSE)
  all_matrix <- lapply(strsplit(chars, "")[[1]], convert_chr_to_matrix, size = LED_height)
  all_com_matrix <- do.call(cbind, all_matrix)

  if (is.null(LED_width)) {
    if (nchar(chars) > 5) {
      LED_width <- 5 * LED_height
    } else {
      LED_width <- nchar(chars) * LED_height
    }
  }
  if (LED_width > ncol(all_com_matrix)) {
    all_com_matrix <- cbind(
      all_com_matrix,
      matrix(0,
        nrow = nrow(all_com_matrix),
        ncol = LED_width - ncol(all_com_matrix)
      )
    )
  }
  width1 <- ncol(all_com_matrix)

  pps <- round(speed / fps)

  pls <- list()
  for (i in seq_len(width1 / pps)) {
    new_all_com_matrix <- cbind(
      all_com_matrix[, (i * pps):width1],
      all_com_matrix[, 1:(i * pps) - 1]
    )
    pls[[i]] <- plot.chr_mat(new_all_com_matrix[, 1:LED_width])
  }

  if (is.null(save_file)) save_file <- file.path(tempdir(), "temp_LED")

  gifski::save_gif(
    {
      for (i in pls) {
        print(i)
      }
    },
    gif_file = paste0(save_file, ".gif"),
    delay = 1 / fps,
    height = 32 * image_scale,
    width = LED_width * image_scale
  )
  magick::image_read(path = paste0(save_file, ".gif"))
}


### 求下一个状态时格子周围值的和：
life_neighbor <- function(m, x, y, size) {
  # m为当前状态的矩阵；x和y为坐标；size为矩阵大小
  fun.sum <- 0
  for (i in c(x - 1, x, x + 1)) { # 依次遍历一个格子周围3x3的邻居格子
    for (j in c(y - 1, y, y + 1)) {
      # 如果格子在角落或者边，则邻居的值直接为0
      if (i > 0 & i <= size & j > 0 & j <= size) fun.sum <- fun.sum + m[i, j] # 把9个格子先求和
    }
  }
  fun.sum <- fun.sum - m[x, y] # 减去中间格子的值，即为周围8个值的和
}

#' Life Game Simulation
#'
#' @param save_file gif filename
#' @param time how many times the life game continue.
#' @param size size of the world
#' @param fps fps, 0.75
#' @param ... add
#' @param colors c("green4", "black")
#'
#' @references <https://zhuanlan.zhihu.com/p/136727731>
#' @return a gif file
#' @export
#' @examples
#' \donttest{
#' if (interactive()) life_game()
#' }
life_game <- function(save_file = NULL, size = 20, time = 20,
                      fps = 0.75, colors = c("black", "green4"), ...) {
  # Game of Life
  ### 构造初始状态：
  # 矩阵的行和列数
  d <- round(runif(size * size, 0, 0.6)) # 最大值低一些，保证初始有值的少一些。
  start <- matrix(data = d, ncol = size, nrow = size)

  ### 设置运行次数
  time <- time
  life <- list()
  life[[1]] <- start
  for (k in 2:time) { # k = 3
    life.next <- matrix(data = 0, ncol = size, nrow = size)
    for (i in 1:size) {
      for (j in 1:size) {
        fun.sum <- life_neighbor(life[[k - 1]], i, j, size)

        # 判断下个状态时当前位置是否有值存在。
        # 孤单死亡：如果细胞的邻居小于等于1个，则该细胞在下一次状态将死亡；
        # 拥挤死亡：如果细胞的邻居在4个及以上，则该细胞在下一次状态将死亡；
        # 稳定：如果细胞的邻居为2个或3个，则下一次状态为稳定存活；
        # 复活：如果某位置原无细胞存活，而该位置的邻居为2个或3个，则该位置将复活一个细胞

        life.next[i, j] <- ifelse(fun.sum == 2 | fun.sum == 3, 1, 0)
      }
    }
    life[[k]] <- life.next
  }
  pls <- lapply(life, plot.chr_mat, colors = colors)

  if (is.null(save_file)) save_file <- file.path(tempdir(), "temp_life")

  gifski::save_gif(
    {
      for (i in pls) {
        print(i)
      }
    },
    gif_file = paste0(save_file, ".gif"),
    delay = 1 / fps,
    width = 400,
    height = 400
  )

  magick::image_read(path = paste0(save_file, ".gif"))
}

#' Plot clock
#'
#' @param x time, defalut: format(Sys.time(), "%H:%M"), e.g. 12:30
#' @param rotate_text rotate_text, FALSE
#' @param text_color text_color, "black"
#' @param bg_color bg_color, "white"
#' @param pointer_color pointer_color, "black"
#' @param time_label time_label, default: as.roman(1:12)
#'
#' @return ggplot
#' @export
#' @references <https://allancameron.github.io/geomtextpath/>
#' @examples
#' clock()
clock <- function(x = format(Sys.time(), "%H:%M"),
                  time_label = as.roman(1:12),
                  rotate_text = FALSE, text_color = "black",
                  bg_color = "white", pointer_color = "black") {
  V1 <- V3 <- V2 <- V4 <- label <- NULL
  if (is.character(x)) x <- as.numeric(strsplit(x, ":")[[1]])
  hours <- c(rep(x[1] %% 12 + tail(x, 1) / 60, 2), 0, 3.5)
  minutes <- c(rep(tail(x, 1) / 5, 2), 0, 5)

  p <- ggplot(as.data.frame(rbind(hours, minutes)), aes(x = V1, y = V3)) +
    geom_rect(xmin = 0, xmax = 12, ymin = 0, ymax = 7.4, fill = bg_color) +
    geom_segment(x = 0, xend = 12, y = 7.4, yend = 7.4, color = "black") +
    geom_segment(
      data = as.data.frame(rbind(hours, minutes)),
      aes(x = V1, y = V3, xend = V2, yend = V4),
      size = c(3, 2), lineend = "round", color = pointer_color
    ) +
    geom_point(x = 0, y = 0, size = 6, color = text_color) +
    scale_x_continuous(limits = c(0, 12), breaks = 1:12) +
    scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) +
    theme_void() +
    coord_polar()
  if (rotate_text) {
    lib_ps("geomtextpath", library = FALSE)
    p <- p +
      geomtextpath::geom_textpath(
        data = data.frame(
          x = 1:12,
          label = as.character(time_label)
        ),
        aes(x = x, y = 6, label = label), size = 10, color = text_color, fontface = "bold"
      )
  } else {
    p <- p +
      geom_text(
        data = data.frame(x = 1:12, label = as.character(time_label)),
        aes(x = x, y = 6, label = label), size = 10, color = text_color, fontface = "bold"
      )
  }
  return(p)
}

#' Plot a github style calendar
#'
#' @param usr github username
#' @param color color, NULL
#' @param save_file save_file, NULL
#' @param ... add
#' @return a svg file
#' @export
show_github_calendar <- function(usr = "asa12138", color = NULL, save_file = NULL, ...) {
  if (is.null(save_file)) {
    save_file <- file.path(tempdir(), paste0(usr, "_calendar.svg"))
  } else {
    if (!grepl(".svg$", save_file)) save_file <- paste0(save_file, ".svg")
  }
  if (!is.null(color)) {
    grDevices::col2rgb(color) %>%
      t() %>%
      grDevices::rgb(., maxColorValue = 255) -> color
  }
  pcutils::download2(paste("https://ghchart.rshah.org", gsub("#", "", color), usr, sep = "/"),
    save_file,
    force = TRUE
  )
  lib_ps("magick", library = FALSE)
  image <- magick::image_read(save_file, density = 300, ...)
  plot(image[1])
}
