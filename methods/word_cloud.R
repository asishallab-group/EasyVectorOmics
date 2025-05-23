# === Load libraries ===
library(dplyr)
library(readr)
library(tm)
library(wordcloud)
library(RColorBrewer)

# === Load outlier gene table ===
distances_with_rdi <- read.table("material/distances_df.tsv", header = TRUE, sep = "\t")

# === Load functional annotations file ===
annotations <- read.table("material/prot_scriber_results.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# === Extract outlier gene IDs ===
outlier_ids <- distances_with_rdi %>%
  filter(is_outlier == TRUE) %>%
  pull(gene_id)

# === Filter annotations for outlier genes ===
outlier_annotations <- annotations %>%
  filter(`Annotee_Identifier` %in% outlier_ids)

# === Build text corpus from descriptions ===
corpus <- Corpus(VectorSource(outlier_annotations$`Human_Readable_Description`))

corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, stripWhitespace)

# === Create Term-Document Matrix ===
tdm <- TermDocumentMatrix(corpus)
tdm_matrix <- as.matrix(tdm)
print(dim(tdm_matrix))
print(head(tdm_matrix))

# === Filter out overly frequent terms (appear in >80% of descriptions) ===
num_docs <- ncol(tdm_matrix)
print(num_docs*0.8)
term_doc_freq <- rowSums(tdm_matrix > 0)
filtered_terms <- term_doc_freq < (num_docs * 0.8)
filtered_matrix <- tdm_matrix[filtered_terms, ]

# === Compute word frequencies ===
word_freq <- sort(rowSums(filtered_matrix), decreasing = TRUE)
word_freq_df <- data.frame(
  word = names(word_freq),
  frequency = as.integer(word_freq)
)
top_20 <- head(word_freq_df, 20)

print(top_20)

write.table(word_freq_df, file = "results/word_freq.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# === Generate word cloud ===
wordcloud(words = names(word_freq),
          freq = word_freq,
          scale = c(4, 0.5),
          max.words = 100,
          random.order = FALSE,
          colors = brewer.pal(8, "Dark2"))

# === Save word cloud to file ===
png("results/outlier_wordcloud_centroid.png", width = 800, height = 600)
wordcloud(words = names(word_freq),
          freq = word_freq,
          scale = c(4, 0.5),
          max.words = 100,
          random.order = FALSE,
          colors = brewer.pal(8, "Dark2"))
dev.off()


# === [NEW] Prepare annotations for all genes ===
# Create corpus for all genes with available annotation
all_annotations <- annotations %>%
  filter(!is.na(Human_Readable_Description))

# === [NEW] Build gene x word presence/absence matrix ===
build_term_matrix <- function(descriptions, gene_ids) {
  corpus <- Corpus(VectorSource(descriptions))
  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, removeNumbers)
  corpus <- tm_map(corpus, stripWhitespace)
  
  tdm <- TermDocumentMatrix(corpus, control = list(wordLengths=c(3, Inf)))  # solo palabras >=3 letras
  tdm_matrix <- as.matrix(tdm)
  colnames(tdm_matrix) <- gene_ids
  return(tdm_matrix)
}

# === [NEW] Generate term matrices ===
all_tdm <- build_term_matrix(all_annotations$Human_Readable_Description, all_annotations$Annotee_Identifier)

# Identify columns corresponding to outlier and non-outlier genes
outlier_cols <- colnames(all_tdm) %in% outlier_ids
non_outlier_cols <- !outlier_cols

# === [NEW] Apply Fisher's exact test to each word ===
fisher_results <- apply(all_tdm, 1, function(word_row) {
  # Contar ocurrencias
  a <- sum(word_row[outlier_cols] > 0)          # outliers with the word
  b <- sum(word_row[outlier_cols] == 0)         # outliers without the word
  c <- sum(word_row[non_outlier_cols] > 0)      # non-outliers with the word
  d <- sum(word_row[non_outlier_cols] == 0)     # non-outliers without the word
  
  contingency <- matrix(c(a, b, c, d), nrow = 2)
  
  # Return NA if contingency table is invalid
  if (any(contingency < 0) || sum(contingency) == 0) {
    return(c(p_value = NA, odds_ratio = NA))
  }
  
  test <- fisher.test(contingency)
  return(c(p_value = test$p.value, odds_ratio = test$estimate))
})

# === [NEW] Convert to data frame and sort by significance ===
fisher_df <- as.data.frame(t(fisher_results))
fisher_df$word <- rownames(fisher_df)
fisher_df <- fisher_df %>%
  arrange(p_value)

# === Save results to TSV file  ===
write.table(fisher_df, file = "results/fisher_keywords.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# === Print top 20 significant terms (p < 0.05) ===
print(head(fisher_df[fisher_df$p_value < 0.05, ], 20))
