from typing import Literal, get_args, get_origin
from pydantic import BaseModel, ConfigDict, Field
from models import RastrModel
import astra

def extract_item_class(annotation):
    """ Находит класс Pydantic-модели внутри сложных типов (Union, List, Optional) """
    # Если это Union (например, list[...] | None), перебираем его аргументы
    if get_origin(annotation) is str or get_origin(annotation).__name__ in ('UnionType', 'Union'):
        for arg in get_args(annotation):
            item_cls = extract_item_class(arg)
            if item_cls:
                return item_cls
    # Если это list[GouItem], достаем сам GouItem
    if get_origin(annotation) is list:
        args = get_args(annotation)
        if args and issubclass(args[0], BaseModel):
            return args[0] 
    return None

def export() :
    rastr = astra.Rastr()
    rastr.load(astra.LoadCode.REPL, "e:/downloads/тест2/тест2/Кубанское РДУ РМ с РИ зима max 2026.os", "d:/documents/rastrwin3\shablon\poisk.os")
    rastr_model = RastrModel.model_construct()
    export_model = {}
    rastr.enum_as_int(False)

    for table_name, table_meta in RastrModel.model_fields.items():
        if not table_name in rastr.tables():
            continue

        rastr_table = rastr.table(table_name)
        item_class = extract_item_class(table_meta.annotation)
        container = []

        for index in range(0, rastr_table.size):
            record = {}
            for field_name, field_meta in item_class.model_fields.items():
                field_name_alias = field_name if not field_meta.alias else field_meta.alias
                rastr_field = rastr_table.column(field_name_alias)
                if rastr_field.type == astra.PropType.ENUM:
                    enum_items = rastr_field.property(astra.FieldProperties.NAMEREF).split("|")
                    string_value = enum_items[rastr_field.z(index)]
                    enum_class = field_meta.annotation
                    enum_instance = enum_class(string_value)
                    record[field_name] = enum_instance
                elif rastr_field.type == astra.PropType.ENPIC:
                    enum_class = field_meta.annotation
                    enum_list = list(enum_class)
                    string_value = enum_list[min(rastr_field.z(index), len(enum_list) - 1)].value
                    enum_instance = enum_class(string_value)
                    record[field_name] = enum_instance
                else:
                    record[field_name] = rastr_field.zn(index)
            container.append(record)
        export_model[table_name] = container
    export_model["version"] = "1.0.0"

    export_rastr_model = RastrModel.model_validate(export_model)

    with open('data.json', 'w', encoding='utf-8') as f:
        f.write(export_rastr_model.model_dump_json(indent=4, by_alias=True))

export()


