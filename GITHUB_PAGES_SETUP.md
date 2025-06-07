# GitHub Pages Setup for fmristore

## Repository Settings Required

To get your pkgdown documentation site working on GitHub Pages, you need to configure the following settings in your GitHub repository:

### 1. Go to Repository Settings
- Navigate to your repository: https://github.com/bbuchsbaum/fmristore
- Click on "Settings" tab (at the top of the repository)

### 2. Configure Pages Settings
- Scroll down to "Pages" in the left sidebar
- Under "Source", select **"GitHub Actions"** (NOT "Deploy from a branch")
- This enables the modern GitHub Pages deployment method

### 3. Verify Workflow Permissions
- Go to "Actions" in the left sidebar
- Click on "General"
- Under "Workflow permissions", ensure **"Read and write permissions"** is selected
- Check "Allow GitHub Actions to create and approve pull requests"

### 4. Environment Settings (Optional but Recommended)
- Go to "Environments" in the left sidebar
- If "github-pages" environment doesn't exist, it will be created automatically
- You can optionally add deployment protection rules here

## After Configuration

1. **Push the updated workflow**: The updated `.github/workflows/pkgdown.yaml` uses the modern GitHub Pages deployment method
2. **Trigger the workflow**: Push any commit to main branch or manually trigger via "Actions" tab
3. **Check deployment**: Go to "Actions" tab to see if the workflow runs successfully
4. **Access your site**: Once deployed, your site will be available at: https://bbuchsbaum.github.io/fmristore/

## Troubleshooting

If the site still doesn't appear:

1. **Check Actions tab** for any failed workflow runs
2. **Verify Pages settings** are set to "GitHub Actions" source
3. **Check repository visibility** - the repo needs to be public for free GitHub Pages
4. **Wait a few minutes** - sometimes there's a delay in DNS propagation

## Current Status

- ✅ Updated workflow to use modern GitHub Pages deployment
- ✅ Added proper permissions for Pages deployment  
- ✅ Configured environment settings for Pages
- ⏳ Need to configure repository settings (see steps above)

Your site URL will be: https://bbuchsbaum.github.io/fmristore/ 