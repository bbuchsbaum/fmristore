# styler configuration for tidyverse style
# This file configures how styler formats code in this package

styler::tidyverse_style(
  scope = I(c("line_breaks", "spaces", "tokens", "indention")),
  strict = FALSE,  # Allow some flexibility for scientific code
  indent_by = 2,   # Standard tidyverse indentation
  start_comments_with_one_space = TRUE,
  reindention = tidyverse_reindention(),
  math_token_spacing = tidyverse_math_token_spacing()
) 