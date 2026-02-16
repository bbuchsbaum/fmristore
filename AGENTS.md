# Repository Guidelines

## Project Structure & Module Organization
- Core package code lives in `R/` (S4 classes, generics, helpers).
- Tests are under `tests/testthat/` with files named `test-*.R`.
- Long-form docs and examples live in `vignettes/`; data and examples used at runtime go in `inst/extdata/` and `inst/examples/`.
- Benchmarks and exploratory scripts belong in `benchmarks/` or `deprecated/`, not in `R/`.

## Build, Test, and Development Commands
- Install dev dependencies: `devtools::install_dev_deps()`.
- Run tests: `devtools::test()`; for a single file, e.g. `devtools::test_file("tests/testthat/test-h5_utils_errors.R")`.
- Full package check: `devtools::check()` before any non-trivial PR.
- Build docs/site: `pkgdown::build_site()` for local documentation preview.
- Run local hooks: `pre-commit run --all-files` to match CI.

## Coding Style & Naming Conventions
- Follow the tidyverse style guide (see `DEVELOPMENT.md`); use 2-space indentation and a soft 80-character line width.
- R objects and functions: `snake_case`; S4 classes: `UpperCamelCase` (e.g., `H5NeuroVol`); tests: `test-<topic>.R`.
- Use roxygen2 comments (`#'`) for all exported functions and classes.
- Before committing, run `styler::style_pkg()` and `lintr::lint_package()`.

## Testing Guidelines
- Tests use `testthat` (edition 3) in `tests/testthat/`.
- Prefer small, focused tests tied to specific files or behaviours.
- Aim for ≥80% coverage; add tests alongside any new feature or bug fix.
- When touching performance-critical paths, add or update cases under `benchmarks/` if relevant.

## Commit & Pull Request Guidelines
- Write concise, present-tense commit messages (e.g., “Add H5 corruption checks”, “Fix cluster summary output”).
- Group related changes into a single commit or small series; avoid mixing refactors with functional changes.
- PRs should describe motivation, key changes, testing performed, and link any relevant issues.
- Ensure CI (R CMD check, style, lint, coverage) passes before requesting review.

## Agent-Specific Notes
- Automated tools and LLM-based agents should make minimal, well-scoped changes, respect this guide and `DEVELOPMENT.md`, and avoid adding new dependencies without prior discussion.
- Prefer updating existing patterns over introducing new styles or architectures.


## Issue Tracking with Beads

This project uses **beads** (`bd`) for git-backed issue tracking. See https://github.com/steveyegge/beads

### Essential Commands

| Command | Purpose |
|---------|---------|
| `bd ready` | List tasks without blockers (your next work) |
| `bd create "title" -p 1` | Create task (P0=critical, P1=high, P2=medium, P3=low) |
| `bd show <id>` | View issue details and history |
| `bd update <id> --status in_progress` | Mark task as in progress |
| `bd close <id> --reason "text"` | Close completed task |
| `bd dep add <child> <parent>` | Add dependency |
| `bd list --json` | List all open issues |
| `bd sync` | Force sync to git |

### Critical Rules for Agents

1. **NEVER use `bd edit`** - it opens an interactive editor. Use flag-based updates:
   ```bash
   bd update <id> --description "new description"
   bd update <id> --title "new title"
   ```

2. **Always use `--json` flag** for programmatic access

3. **Run `bd sync` after changes** to ensure immediate git sync

### Finding Work

```bash
bd ready --json          # Tasks without blockers
bd list --status open    # All open tasks
bd stale --days 7        # Neglected tasks
```

## Landing the Plane (Session Completion)

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd sync
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
