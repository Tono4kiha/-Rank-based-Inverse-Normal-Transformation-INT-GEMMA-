#!/user/bin/Rscript

# ============================================================
# 表型 Rank-based Inverse Normal Transformation (INT)
# 输出 GEMMA 可直接使用的 phenotype 文件
#
# 用法:
#   Rscript pheno_int.R <文件>
#   Rscript pheno_int.R <文件> trait1 trait2
#   Rscript pheno_int.R <文件> yes
#   Rscript pheno_int.R <文件> trait1 trait2 yes
#
# 思路:
#   rank -> percentile -> qnorm
#   使用 Blom transformation:
#
#   z = qnorm((r - 3/8) / (N + 1/4))
#
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {

  cat("用法:\n")
  cat("  Rscript pheno_int.R pheno.txt\n")
  cat("  Rscript pheno_int.R pheno.txt trait1 trait2\n")
  cat("  Rscript pheno_int.R pheno.txt yes\n")
  cat("  Rscript pheno_int.R pheno.txt trait1 trait2 yes\n")

  quit(status = 1)
}

input_file <- args[1]
rest <- args[-1]

# ------------------------------------------------------------
# 判断最后一个参数是否为 yes/no
# ------------------------------------------------------------

plot_flag <- "no"

if (length(rest) > 0 &&
    tail(rest, 1) %in% c("yes", "no")) {

  plot_flag <- tail(rest, 1)
  rest <- head(rest, -1)
}
user_pheno_names <- rest
# ------------------------------------------------------------
# 自动判断是否有表头
# ------------------------------------------------------------

first_line <- readLines(input_file, n = 1)

fields <- strsplit(
  trimws(first_line),
  "\\s+"
)[[1]]

has_header <- any(
  is.na(
    suppressWarnings(
      as.numeric(fields[-1])
    )
  )
)

# ------------------------------------------------------------
# 读入数据
# ------------------------------------------------------------

# -------- 读取数据并设置列名 --------
df <- read.table(
  input_file,
  header = has_header,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

n_pheno_in_file <- ncol(df) - 1

if (has_header) {

  # 有表头：保留文件中的表型列名，仅统一第一列为 id
  colnames(df)[1] <- "id"

} else {

  # 无表头但用户输入了表型名：
  # 按用户输入顺序给第2列及之后的表型列命名
  if (length(user_pheno_names) > 0) {

    if (length(user_pheno_names) != n_pheno_in_file) {
      stop(
        sprintf(
          paste0(
            "检测到输入文件没有表头。\n",
            "文件中共有 %d 个表型列，但命令中提供了 %d 个表型名称：%s\n",
            "无表头模式下，输入的表型名称数量必须与文件中的表型列数完全一致。"
          ),
          n_pheno_in_file,
          length(user_pheno_names),
          paste(user_pheno_names, collapse = ", ")
        )
      )
    }

    colnames(df) <- c("id", user_pheno_names)

  } else {

    # 无表头且未提供名称：使用默认名称
    colnames(df) <- c(
      "id",
      paste0("trait", seq_len(n_pheno_in_file))
    )
  }
}
# ------------------------------------------------------------
# 确定表型列
# ------------------------------------------------------------
if (has_header) {

  if (length(user_pheno_names) == 0) {

    # 有表头且未指定表型：转换除 id 外的全部表型
    pheno_cols <- setdiff(colnames(df), "id")

  } else {

    # 有表头且指定了表型：按照列名选择
    pheno_cols <- user_pheno_names

    miss <- setdiff(pheno_cols, colnames(df))

    if (length(miss) > 0) {
      stop(
        sprintf(
          "列名不存在: %s\n可选列: %s",
          paste(miss, collapse = ", "),
          paste(setdiff(colnames(df), "id"), collapse = ", ")
        )
      )
    }
  }

} else {

  # 无表头时，列名已经在上一部分完成设置
  # 默认转换文件中的全部表型列
  pheno_cols <- setdiff(colnames(df), "id")
}

cat("\n========== 输入信息 ==========\n")
cat("文件       :", input_file, "\n")
cat("表头       :", ifelse(has_header, "有", "无"), "\n")
cat("数据列名   :", paste(colnames(df), collapse = ", "), "\n")
cat("转换表型   :", paste(pheno_cols, collapse = ", "), "\n")
cat("输出检查图 :", plot_flag, "\n")
# ------------------------------------------------------------
# 缺失值处理
# ------------------------------------------------------------

missing_vals <- c(
  -9,
  "-9",
  "NA",
  "na",
  "NaN",
  "."
)

for (col in pheno_cols) {

  x <- df[[col]]

  if (is.character(x)) {
    x <- suppressWarnings(as.numeric(x))
  }

  x[x %in% missing_vals] <- NA

  df[[col]] <- x
}

# ------------------------------------------------------------
# 表型QC
# ------------------------------------------------------------

cat("\n========== 表型质量检查 ==========\n")

for (col in pheno_cols) {

  x <- df[[col]]

  n_valid <- sum(!is.na(x))
  n_miss <- sum(is.na(x))
  n_unique <- length(unique(x[!is.na(x)]))

  cat(
    sprintf(
      "%-20s N=%d  Missing=%d  Unique=%d\n",
      col,
      n_valid,
      n_miss,
      n_unique
    )
  )

  if (n_valid < 30) {

    warning(
      sprintf(
        "%s 有效样本数仅 %d，GWAS统计功效可能不足",
        col,
        n_valid
      )
    )
  }

  if (n_unique == 1) {

    warning(
      sprintf(
        "%s 为常数表型（所有值相同）",
        col
      )
    )
  }

  if (n_unique > 1 &&
      n_unique < 5) {

    warning(
      sprintf(
        "%s 仅有 %d 个不同取值，请确认是否适合连续性状GWAS",
        col,
        n_unique
      )
    )
  }
}

# ------------------------------------------------------------
# Rank-based INT
# Blom transformation
# ------------------------------------------------------------

rank.based.INT <- function(
  x,
  c = 3/8,
  method = "average"
) {

  r <- rank(
    x,
    ties.method = method
  )

  r[is.na(x)] <- NA

  N <- sum(!is.na(x))

  qnorm(
    (r - c) /
      (N - 2 * c + 1)
  )
}

# ------------------------------------------------------------
# INT转换并输出GEMMA文件
# ------------------------------------------------------------

cat("\n========== 转换输出 ==========\n")

for (col in pheno_cols) {

  x <- df[[col]]

  n_valid <- sum(!is.na(x))
  n_unique <- length(unique(x[!is.na(x)]))

  if (n_valid < 3) {

    warning(
      sprintf(
        "%s 有效样本数不足，跳过",
        col
      )
    )

    next
  }

  if (n_unique == 1) {

    warning(
      sprintf(
        "%s 为常数表型，跳过",
        col
      )
    )

    next
  }

  int_col <- paste0(
    col,
    "_INT"
  )

  df[[int_col]] <- rank.based.INT(x)

  out_file <- sprintf(
    "pheno_%s_GEMMA.txt",
    col
  )

  write.table(

    data.frame(
      FID = df$id,
      IID = df$id,
      Pheno = df[[int_col]]
    ),

    file = out_file,

    col.names = FALSE,
    row.names = FALSE,

    quote = FALSE,
    sep = "\t",

    na = "NA"
  )

  cat(
    sprintf(
      "  ✓ %-20s -> %s\n",
      col,
      out_file
    )
  )
}

# ------------------------------------------------------------
# 可选检查图
# ------------------------------------------------------------

if (plot_flag == "yes") {

  out_jpg <- "INT_check.jpg"

  jpeg(
    out_jpg,
    width = 1200,
    height = 800
  )

  par(
    mfrow = c(
      length(pheno_cols),
      2
    )
  )

  for (col in pheno_cols) {

    int_col <- paste0(
      col,
      "_INT"
    )

    if (!int_col %in% colnames(df)) {
      next
    }

    hist(
      df[[int_col]],
      breaks = 30,
      main = paste(col, "(INT)"),
      col = "lightgreen",
      xlab = col
    )

    qqnorm(
      df[[int_col]],
      main = paste("QQ:", col)
    )

    qqline(
      df[[int_col]],
      col = "red"
    )
  }

  dev.off()

  cat(
    "\n检查图: ",
    out_jpg,
    "\n"
  )
}

cat("\n✅ 完成\n")
``
