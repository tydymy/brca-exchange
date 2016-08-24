# -*- coding: utf-8 -*-
# Generated by Django 1.9.2 on 2016-02-10 16:52
from __future__ import unicode_literals

import csv
import os.path

from django.conf import settings
from django.db import migrations

from data.models import Variant


def load_from_csv(apps, schema_editor):
    file_path = os.path.join(settings.BASE_DIR, 'data', 'resources', 'written.tsv')
    with open(file_path) as tsv_file:
        reader = csv.reader(tsv_file, dialect="excel-tab")
        header = reader.next()

        for row in reader:

            # split Source column into booleans
            row_dict = dict(zip(header, row))
            for source in row_dict['Source'].split(','):
                row_dict['Variant_in_' + source] = True
            Variant.objects.create_variant(row_dict)


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0002_search_index'),
    ]

    operations = [
        migrations.RunSQL("DELETE FROM variant;"),
        migrations.RunPython(load_from_csv),
    ]
