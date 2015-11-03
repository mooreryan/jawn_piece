module Const
  WINDOW_SIZE = 300
  WINDOW_STEP = 25

  SILVA_LEN = 50000

  DATABASE = File.join "test_files", "database.fa"

  ASSETS_FOLDER = "assets"

  DE_DIST = File.join ASSETS_FOLDER, "database_DE_dist.txt"
  DE_DIST_PERCENTILES =
    File.join ASSETS_FOLDER, "database_DE_dist_percentiles.txt"

  NUM_PERCENTILES = 10

  RSCRIPT = File.join ASSETS_FOLDER, "pintail_plot.r"

end

Ryan.try_mkdir Const::ASSETS_FOLDER
