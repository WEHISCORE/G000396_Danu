# Copy large files not included in git repo to OneDrive.
# Peter Hickey
# 2023-12-05

# Project specific variables ---------------------------------------------------

REMOTE="WEHIOneDrive"
PROJECT="G000396_Danu"
PROJECT_ROOT="/stornext/Projects/score/Analyses/${PROJECT}"
LOCAL_PATH="${PROJECT_ROOT}"
REMOTE_PATH="SCORE/${PROJECT}"

# Create destination directory (it it doesn't exist already) -------------------

rclone mkdir ${REMOTE}:${REMOTE_PATH}

# Copy directories and their contents from local to remote ---------------------

DIRS=( DEGs Glimma heatmaps scPipe )
for DIR in "${DIRS[@]}"
do
  echo "Copying ${DIR}"
  rclone copy --progress ${LOCAL_PATH}/output/${DIR} ${REMOTE}:${REMOTE_PATH}/output/${DIR}
  echo ""
done
