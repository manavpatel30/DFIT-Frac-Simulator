# Publish the Project on GitHub

## Option 1: upload in a browser

1. Sign in to GitHub and create a new repository.
2. Leave **Add a README**, **Add .gitignore**, and **Choose a license** unchecked during repository creation because this project already contains a README and `.gitignore`.
3. Extract the provided ZIP file.
4. On the empty repository page, choose **uploading an existing file**.
5. Drag all contents from inside the extracted `dfit-fracture-simulator` folder into the upload area. Include the hidden `.gitignore` file.
6. Enter a commit message such as `Initial DFIT simulator release` and commit the files.

Choose a license separately before allowing others to reuse or modify the project.

## Option 2: use Git from a terminal

Create an empty GitHub repository first. Then run the following commands inside the extracted project folder, replacing the placeholder URL with the URL of your repository:

```bash
git init
git add .
git commit -m "Initial DFIT simulator release"
git branch -M main
git remote add origin https://github.com/USERNAME/REPOSITORY.git
git push -u origin main
```

## Suggested repository description

> Physics-based Python simulator and interactive dashboard for a simplified sleeve-operated DFIT, including pressurization, PKN fracture propagation, Carter leakoff, friction losses, and post-shut-in pressure decline.

## Suggested repository topics

`dfit`, `hydraulic-fracturing`, `geomechanics`, `fracture-mechanics`, `petroleum-engineering`, `pkn-model`, `python`, `simulation`

