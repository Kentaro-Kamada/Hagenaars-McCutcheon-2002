library(tidyverse)
library(kamaken)
library(gt)
set.seed(123)

# データの読み込み --------------------------------------------------------------------------------------


'Score	B	C	n
0	0	0	61
0	0	1	24
0	1	0	9
0	1	1	6
1	0	0	92
1	0	1	50
1	1	0	28
1	1	1	17
2	0	0	107
2	0	1	60
2	1	0	30
2	1	1	37
3	0	0	75
3	0	1	44
3	1	0	32
3	1	1	71
4	0	0	56
4	0	1	30
4	1	0	44
4	1	1	102
5	0	0	28
5	0	1	17
5	1	0	20
5	1	1	112
6	0	0	5
6	0	1	2
6	1	0	9
6	1	1	73' %>% 
  write_file('.tmp')

data <- read_tsv('.tmp')

file.remove('.tmp')

# 個票データの形式にする
df <- 
  array(data$n, dim = c(2, 2, 7), 
        dimnames = list(C = unique(data$C), 
                        B = unique(data$B), 
                        Score = unique(data$Score))) %>% 
  epitools::expand.table() %>% 
  as_tibble() %>% 
  relocate(B, C, Score) %>%  
  mutate(across(c(B, C), .fns = as.numeric)) %>% 
  mutate(Score = as.character(Score) %>% parse_double())
  # mutate(Score = as.numeric(Score))

# データ書き出し
write_delim(df, 'geometry.dat', delim = ' ', col_names = F, eol = '\r\n')
read_file('geometry.dat') %>% 
  str_remove('\r\n$') %>% 
  write_file('geometry.dat')



# lemの実行 --------------------------------------------------------------------------------------


lem <- function(lat, man, con, dim, lab, mod, rec, des, dum, dat, ite, see, path) {
  # inpとoutのディレクトリ作成
  if(!dir.exists(path)) {
    dir.create(path)
    dir.create(str_c(path, '/inp'))
    dir.create(str_c(path, '/out'))
  }
  
  # ファイルの削除
  file.remove(list.files(path, full.names = T, recursive = T))
  
  # スクリプト書き出し
  walk(see, ~{write_lines(c(lat, man, con, dim, lab, mod, rec, des, dum, dat, ite, .), 
                          file = str_c(path, 'inp/', str_remove(., ' '), '.inp'))})
  
  
  # inputファイルとoutputファイルのパス指定
  inp_path <- list.files(str_c(path, 'inp'), full.names = T)
  out_path <- str_replace_all(inp_path, 'inp', 'out')
  
  # lem実行
  walk(str_c('lem95', inp_path, out_path, sep = ' '),
       ~{system(., wait = F)
         Sys.sleep(0.4)})
  
  # 最大対数尤度の読み出し
  max_loglik <-
    tibble(
      filepath = out_path,
      loglik = map_chr(filepath, read_file)
    ) %>% 
    mutate(
      seed = str_extract(loglik, '(?<=Seed random values   = ).+') %>% 
        parse_double(),
      loglik = str_extract(loglik, '(?<=Log-likelihood       = ).+') %>% 
        parse_double()
    ) %>% 
    arrange(desc(loglik)) %>% 
    mutate(rank = dense_rank(desc(loglik)))
  
  return(max_loglik)
}

lem_max_loglik <- function(out_path) {
  
  # 最大対数尤度の読み出し
  tibble(
    filepath = list.files(out_path, full.names = T),
    loglik = map_chr(filepath, read_file)
  ) %>% 
    mutate(
      seed = str_extract(loglik, '(?<=Seed random values   = ).+') %>% 
        parse_double(),
      loglik = str_extract(loglik, '(?<=Log-likelihood       = ).+') %>% 
        parse_double()
    ) %>% 
    arrange(desc(loglik)) %>% 
    mutate(rank = dense_rank(desc(loglik)))
}



# tab7 --------------------------------------------------------------------------------------
# 使用データに注意


model1 <- 
  lem(lat = 'lat 1',
      man = 'man 2',
      con = 'con 1',
      dim = 'dim 2 2 2',
      lab = 'lab X A B x', 
      mod = 'mod X|x {X, cov(x, 1, X, c, -1)} 
                 A|X
                 B|X',
      des = 'des [0 1]',
      dum = 'dum 1 1 1',
      rec = 'rec 1241', 
      ite = 'ite 10000',
      dat = 'dat geometry.dat',
      see = str_c('see ', rdunif(20, 999, 100)),
      path = 'lem/tab7/model1/'
  )

model2 <- 
  lem(lat = 'lat 1',
      man = 'man 2',
      con = 'con 1',
      dim = 'dim 2 2 2',
      lab = 'lab X A B x', 
      mod = 'mod X|x {X, cov(x, 1, X, c, -1)} 
                 A|X eq2
                 B|X eq2',
      des = 'des [0 1
                  1 0 2 0
                  1 0 2 0
      ]',
      dum = 'dum 1 1 1',
      rec = 'rec 1241', 
      ite = 'ite 10000',
      dat = 'dat geometry.dat',
      see = str_c('see ', rdunif(20, 999, 100)),
      path = 'lem/tab7/model2/'
  )

model3 <- 
  lem(lat = 'lat 1',
      man = 'man 3',
      con = '',
      dim = 'dim 2 2 2 7',
      lab = 'lab X A B x', 
      mod = 'mod X|x 
                 A|X
                 B|X',
      des = '',
      dum = 'dum 1 1 1 1',
      rec = 'rec 1241', 
      ite = 'ite 10000',
      dat = 'dat geometry.dat',
      see = str_c('see ', rdunif(20, 999, 100)),
      path = 'lem/tab7/model3/'
  )

model4 <- 
  lem(lat = 'lat 1',
      man = 'man 3',
      con = '',
      dim = 'dim 2 2 2 7',
      lab = 'lab X A B x', 
      mod = 'mod X|x 
                 A|X eq2
                 B|X eq2',
      des = 'des [1 0 2 0
                  1 0 2 0]',
      dum = 'dum 1 1 1 1',
      rec = 'rec 1241', 
      ite = 'ite 10000',
      dat = 'dat geometry.dat',
      see = str_c('see ', rdunif(20, 999, 100)),
      path = 'lem/tab7/model4/'
  )

model1 <- lem_max_loglik('lem/tab7/model1/out/')
model2 <- lem_max_loglik('lem/tab7/model2/out/')
model3 <- lem_max_loglik('lem/tab7/model3/out/')
model4 <- lem_max_loglik('lem/tab7/model4/out/')


# 適合度指標回収

tibble(
  model = str_c('model', 1:4),
  filepath = c(model1$filepath[1], model2$filepath[1], model3$filepath[1], 
               model4$filepath[1]),
  output = map_chr(filepath, read_file)
) %>% 
  transmute(
    model,
    Lsq = str_extract(output, 
                      '(?<=L-squared            = ).+(?=\\(.+\\))') %>% 
      parse_double(),
    Xsq = str_extract(output, 
                      '(?<=X-squared            = ).+(?=\\(.+\\))') %>% 
      parse_double(),
    df = str_extract(output, '(?<=Degrees of freedom   = ).+') %>% 
      parse_double(),
    AIC = str_extract(output, '(?<=AIC\\(L-squared\\)       = ).+') %>% 
      parse_double()
  ) %>% 
  gt() %>% 
  fmt_number(columns = c(Lsq, Xsq, AIC), decimals = 3)
