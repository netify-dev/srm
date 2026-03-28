# generate example datasets for the srm package
#
# these seeds produce the specific values used in all vignette narratives
# and saved .rda files; changing them will break vignette output

# -------------------------------------------------------------------
# classroom: mean girls friendship ratings (hand-designed)
# -------------------------------------------------------------------
set.seed(6886)
n = 12
actors = c("Cady", "Aaron", "Gretchen", "Karen", "Janis",
           "Damian", "Kevin", "Glen", "Shane", "Trang",
           "Regina", "Norbury")

# actor effects: sender generosity (damian warmest, regina stingiest)
a = c(
  0.6,    # cady: kind new girl
  0.2,    # aaron: moderately generous
  -0.3,   # gretchen: catty despite trying
  1.2,    # karen: sweet, naive, rates everyone well
  -0.1,   # janis: selective, prickly to non-friends
  2.0,    # damian: warmest rater in the school
  0.0,    # kevin: neutral
  -0.2,   # glen: barely involved
  -1.0,   # shane: jock, low effort
  -0.4,   # trang: peripheral
  -2.4,   # regina: queen bee, rates strategically low
  0.4     # norbury: supportive teacher
)

# partner effects: receiver popularity (regina most popular, gretchen least)
b = c(
  0.8,    # cady: increasingly popular
  1.3,    # aaron: popular, well-liked
  -1.6,   # gretchen: "none for gretchen wieners"
  0.3,    # karen: liked as a plastic
  -0.4,   # janis: outsider
  -0.3,   # damian: not cool-kid popular
  -0.7,   # kevin: mathlete nerd
  0.1,    # glen: "you go glen coco" mild positive
  -0.4,   # shane: average jock
  -0.4,   # trang: peripheral
  1.5,    # regina: universally feared/admired
  -0.2    # norbury: teacher, not student-popular
)

# dyadic noise with reciprocity
g = matrix(rnorm(n * n, sd = 1.2), n, n)
g = (g + t(g) * 0.4) / 1.4
diag(g) = 0

# hand-set key relationships from the movie
g[1, 5] = 1.5; g[5, 1] = 1.5   # cady-janis: genuine friendship
g[1, 6] = 1.0; g[6, 1] = 1.0   # cady-damian: friendship
g[1, 2] = 1.5; g[2, 1] = 1.2   # cady-aaron: romantic tension
g[1, 11] = 1.0; g[11, 1] = 0.5 # cady-regina: infiltrating / recruiting
g[3, 11] = 2.5; g[11, 3] = -1.0 # gretchen worships regina; regina uses her
g[4, 3] = 0.8; g[3, 4] = 0.8   # karen-gretchen: plastics bond
g[4, 11] = 0.5; g[11, 4] = -0.3 # karen-regina: oblivious plastic
g[5, 11] = -2.0; g[11, 5] = -1.0 # janis-regina: enemies
g[5, 6] = 1.5; g[6, 5] = 1.5   # janis-damian: best friends
g[9, 11] = 1.5; g[11, 9] = 0.3 # shane-regina: romantic
g[7, 1] = 1.5; g[1, 7] = 0.3   # kevin-cady: crush
g[12, 1] = 0.8; g[1, 12] = 0.5 # norbury-cady: mentor

Y = 3.5 + matrix(a, n, n) + matrix(b, n, n, byrow = TRUE) + g
diag(Y) = 0
Y = round(Y, 2)
rownames(Y) = colnames(Y) = actors

classroom = Y

# -------------------------------------------------------------------
# trade_net: simulated bilateral trade intensity (3 periods)
# -------------------------------------------------------------------
set.seed(2025)
countries = c("USA", "CHN", "DEU", "JPN", "GBR", "FRA", "KOR", "IND",
              "BRA", "CAN")
nc = length(countries)
years = c("2015", "2017", "2019")

trade_net = list()
for (yr in years) {
  a = rnorm(nc, sd = 0.9)
  b = rnorm(nc, sd = 0.7)
  g = matrix(rnorm(nc * nc, sd = 1.0), nc, nc)
  g = (g + t(g) * 0.3) / 1.3
  diag(g) = 0
  Y = 2 + matrix(a, nc, nc) + matrix(b, nc, nc, byrow = TRUE) + g
  diag(Y) = 0
  Y = round(Y, 2)
  rownames(Y) = colnames(Y) = countries
  trade_net[[yr]] = Y
}

# -------------------------------------------------------------------
# small_council: GoT houses pursuing small council positions (bipartite)
# -------------------------------------------------------------------
set.seed(6886)

houses = c("Stark", "Lannister", "Targaryen", "Baratheon", "Tyrell",
           "Martell", "Greyjoy", "Arryn", "Tully", "Bolton")
positions = c("Hand", "Coin", "Whispers", "Ships", "War", "Law", "Faith")

nh = length(houses)
np = length(positions)

# actor effects: overall ambition level
a = c(-0.3, 2.8, 2.2, 0.8, 0.5, -0.2, -1.0, -1.5, -1.2, -0.5)

# partner effects: position desirability
b = c(2.0, 1.0, 0.3, -0.5, 0.8, -1.2, -2.0)

# dyadic noise + deliberate relationship-specific patterns
g = matrix(rnorm(nh * np, sd = 0.8), nh, np)
g[1, 1] = 2.0   # stark wants hand (ned)
g[1, 3] = -2.5  # stark avoids whispers (honor)
g[2, 2] = 2.0   # lannister wants coin (gold)
g[2, 5] = 1.5   # lannister wants war (military)
g[3, 1] = 1.5   # targaryen wants hand (rule)
g[3, 5] = 1.8   # targaryen wants war (dragons)
g[5, 2] = 2.5   # tyrell wants coin (wealth)
g[5, 7] = 1.5   # tyrell wants faith (margaery)
g[6, 3] = 2.8   # martell wants whispers (oberyn)
g[7, 4] = 3.5   # greyjoy wants ships (ironborn)
g[7, 1] = -2.0  # greyjoy avoids hand (no interest in court)
g[10, 5] = 2.5  # bolton wants war (ramsay)
g[10, 3] = 1.5  # bolton wants whispers (scheming)
g[8, 6] = 1.5   # arryn wants law (honor)
g[9, 7] = 1.0   # tully wants faith (family values)

Y = 5.0 + matrix(a, nh, np) + matrix(b, nh, np, byrow = TRUE) + g
Y = round(pmax(Y, 0.5), 1)
rownames(Y) = houses
colnames(Y) = positions

small_council = Y

# save
usethis::use_data(classroom, overwrite = TRUE)
usethis::use_data(trade_net, overwrite = TRUE)
usethis::use_data(small_council, overwrite = TRUE)
