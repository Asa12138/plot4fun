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
#' my_wordcloud(taxonomy$Genus)
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
#' @references https://mp.weixin.qq.com/s/W-BYPR3UXL120XWpTmN3rA
give_you_a_rose <- function(color = "red3") {
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
#' @references \code{https://github.com/SherryDong/create_plot_by_R_base}
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
#'
#' @return chr_mat
#' @export
convert_img_to_matrix <- function(image_file, size = 32) {
    lib_ps("magick", library = FALSE)
    # 读取并返回图像
    image <- magick::image_read(image_file)
    # 将图像转换为灰度并调整大小
    image <- magick::image_convert(image, format = "png", depth = 8)
    image <- magick::image_scale(image, paste0(size, "x", size))

    # 提取像素值并转换为01矩阵
    chr_mat <- matrix(as.integer(magick::image_data(image)[1, , ] < 128), nrow = size, byrow = TRUE)
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
#'
#' @return plot
#' @exportS3Method
#' @method plot chr_mat
plot.chr_mat <- function(x, colors = c("grey", "red2"), ...) {
    d <- x
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
#' @param speed 2
#' @param ... add
#' @param LED_width LED_width, 64
#'
#' @return gif file
#' @export
#'
make_LED <- function(chars = "HHH", save_file = NULL,
                     LED_width = 64, speed = 2, colors = c("grey", "red2"), ...) {
    lib_ps("animation", library = FALSE)
    all_matrix <- lapply(strsplit(chars, "")[[1]], convert_chr_to_matrix)
    all_com_matrix <- do.call(cbind, all_matrix)

    pls <- list()
    speed1 <- 2
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
    for (i in seq_len(width1 / speed1)) {
        new_all_com_matrix <- cbind(
            all_com_matrix[, (i * speed1):width1],
            all_com_matrix[, 1:(i * speed1) - 1]
        )
        pls[[i]] <- plot.chr_mat(new_all_com_matrix[, 1:LED_width])
    }

    if (is.null(save_file)) save_file <- file.path(tempdir(), "temp_LED")
    plotgif(pls, file = save_file, speed = speed, height = 32 * 10, width = LED_width * 10)

    magick::image_read(path = paste0(save_file, ".gif"), ...)
}