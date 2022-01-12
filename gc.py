import os

from google.cloud import storage
from google.cloud.storage import constants


def get_paths(client, old_bucket_dir):
    counter = 0
    all_paths = dict()

    for blob in client.list_blobs(old_bucket_dir, prefix='2021-12-13'):
        folder = os.path.dirname(blob.name)
        if 'illumina' not in folder:
            continue
        else:
            folder_spl = folder.split('/')
            if len(folder_spl) == 3:
                key_dir = folder_spl[1]
                if key_dir not in all_paths:
                    all_paths[key_dir] = list()
            elif len(folder_spl) > 3 and '.tar.gz' not in blob.name:
                all_paths[key_dir].append(blob)

            if 'results' not in folder:
                print(folder)

        counter += 1
        if counter == 10:
            break

    return all_paths


def copy_1(source_bucket, destination_bucket):
    counter = 0

    for key in all_paths:
        print(len(all_paths[key]))
        print(all_paths[key])
        for blob in all_paths[key]:

            new_blob_name = '/'.join(blob.name.split('/')[1:])
            if blob.name.split('/')[-1].split('.')[1] in ['annot', 'bam', 'coverage', 'vcf']:
                _ = source_bucket.copy_blob(blob, destination_bucket, new_blob_name)
            # elif '_consensus' in blob.name or '_filtered' in blob.name:
            #    print(blob.name, new_blob_name)
            #    new_blob_name_spl = new_blob_name.split('/')
            #    new_blob_name_spl[2] = f"{new_blob_name_spl[2]}.tar.gz"
            #    new_blob_name = '/'.join(new_blob_name_spl)
            #    _ = source_bucket.copy_blob(blob, destination_bucket, new_blob_name)

            counter += 1
            if counter == 10:
                break

        break

    return


if __name__ == "__main__":
    client = storage.Client()

    old_bucket_dir = 'prj-int-dev-covid19-nf-gls'
    new_bucket_dir = '2022-01-25-upd-bucket'

    all_paths = get_paths(client=client, old_bucket_dir=old_bucket_dir)

    print(all_paths.keys())
    print(len(all_paths))
    for key in all_paths:
        print(key, len(all_paths[key]))

    new_bucket = client.bucket(new_bucket_dir)
    if not new_bucket.exists():
        new_bucket.create()

    source_bucket = client.get_bucket(old_bucket_dir)
    destination_bucket = client.get_bucket(new_bucket_dir)

    copy_1(source_bucket=source_bucket, destination_bucket=destination_bucket)

    destination_bucket.storage_class = constants.ARCHIVE_STORAGE_CLASS

    blobs_names = list()
    for blob in client.list_blobs(new_bucket_dir):
        print(blob.name)
        blob_name = "/".join(blob.name.split('/')[:3])
        if blob_name not in blobs_names:
            blobs_names.append(blob_name)
            # blob.update_storage_class("ARCHIVE")
            metadata = {'Content-Type': 'application/octet-stream', 'Content-Encoding': 'gzip'}
            blob.metadata = metadata
            blob.patch()
