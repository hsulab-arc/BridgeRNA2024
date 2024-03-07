export PROJECT=$(gcloud config get-value project)
export REPO=bridge-analysis
export TAG=latest
export IMAGE_URI=gcr.io/$PROJECT/$REPO:$TAG
docker build . --tag $IMAGE_URI
