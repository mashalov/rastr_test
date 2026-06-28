from ast import Pass
import json
from datamodel_code_generator import InputFileType, generate
from datamodel_code_generator.format import Formatter
from datamodel_code_generator.types import PythonVersion
import re
import astra


mode_from_rastr = True

export_limits = {
    "tables": {
        "exclude" : ["DFWReferenceValues"]
        }
}


def transliterate(text: str) -> str:
    """Переводит кириллицу в латиницу и делает строку валидным ID для Python."""
    cyrillic = "абвгдеёжзийклмнопрстуфхцчшщъыьэюя"
    latin = "abvgdeezzijklmnoprstufhcshsh''eua"
    tr_map = str.maketrans(cyrillic, latin)

    text = text.lower().translate(tr_map)
    # Заменяем пробелы, дефисы и спецсимволы на нижнее подчеркивание
    text = re.sub(r"[^a-z0-9_]", "_", text)
    # Если строка начинается с цифры, добавляем префикс val_
    if text and text[0].isdigit():
        text = f"val_{text}"
    return text or "empty"

enpic_map = {
    ("vetv","sta"): {
        "enum_items" : ["Вкл", "Откл", "Откл в начале", "Откл в конце"],
        "enum_keys" : ["On", "Off", "TripHead", "TripTail"],
        "enum_name" : "BranchStateEnum"
        },
     ("vetv","signP"): {
        "enum_items" : ["От шин", "В шины"],
        "enum_keys" : ["HeadTail", "TailHead"],
        "enum_name" : "BranchActivePowerDirectionEnum"
        },
     ("vetv","signQip"): {
        "enum_items" : ["От шин", "В шины"],
        "enum_keys" : ["HeadTail", "TailHead"],
        "enum_name" : "BranchHeadReactivePowerDirectionEnum"
        },
     ("vetv","signQiq"): {
        "enum_items" : ["От шин", "В шины"],
        "enum_keys" : ["HeadTail", "TailHead"],
        "enum_name" : "BranchTailReactivePowerDirectionEnum"
        }
    }

def generate_model() -> json:
    rastr = astra.Rastr()
    json_schema = {
    "$schema": "http://json-schema.org",
    "title": "RastrModel",
    "type": "object",
    "x-pydantic": {"populate_by_name": True}, 
    "properties": {
        "version": {
                "type": "string",
                "default": "1.0.0",
                "const": "1.0.0",
                "description": "Версия модели"
            }
        },
    "$defs": {},
    "additionalProperties": False,
    "required": ["version"]
 }
    type_map = {
        astra.PropType.INT: "integer", 
        astra.PropType.DOUBLE: "number",
        astra.PropType.STRING: "string",
        astra.PropType.BOOL: "boolean",
        }

    rastr.new_file("d:/documents/rastrwin3\shablon\poisk.os")
    for table in rastr.tables():

        if table.name in export_limits["tables"]["exclude"]:
            continue

        table_keys = table.key.split(",")
        table_data = {
            "type": "object",
            "additionalProperties": False
        }
        fields_data = {}
        required = []
        for field in table.columns():

            if field.name.startswith("_"):
                continue

            add_to_required = field.name in table_keys

            read_only_field = len(field.property(astra.FieldProperties.EXPRESSION)) > 0

            if mode_from_rastr:
                # отдаем вычисляемые поля
                add_to_required = True
            else:
                # не принимаем вычисляемые поля
                if read_only_field:
                    continue

            field_data = {}

            field_description = field.property(astra.FieldProperties.DESCRIPTION)
            field_title = field.property(astra.FieldProperties.TITLE)
            if field_description:
                field_data["description"] = field_description
            if field_title:
                field_data["title"] = field_title

            if field.type in type_map:
                field_data["type"] = type_map[field.type]
            else:
                enum_items = None
                enum_name = None
                enum_keys = None
                match field.type:
                    case astra.PropType.ENUM:
                        enum_items = field.property(astra.FieldProperties.NAMEREF).split("|")
                        enum_name = f"{table.name}{field.name}Enum"
                        enum_keys = [transliterate(x) for x in enum_items]
                    case astra.PropType.ENPIC:
                        enpic_info = enpic_map.get((table.name, field.name), None)
                        if enpic_info:
                            enum_items = enpic_info["enum_items"]
                            enum_name = enpic_info["enum_name"]
                            enum_keys = enpic_info["enum_keys"]
                        else:
                            print(f"Enpic {table.name}.{field.name}")
                    case astra.PropType.SUPERENUM:
                        pass
                        #print(f"Superenum {table.name}.{field.name}")

                if enum_items:
                    json_schema["$defs"][enum_name] = {
                            "type": "string",
                            "title": enum_name,
                            "enum": enum_items,
                            "x-enum-varnames": enum_keys
                        }
                    field_data["$ref"] = f"#/$defs/{enum_name}"

            fields_data[field.name] = field_data
            add_to_required |= field.type == astra.PropType.BOOL
            add_to_required |= len(field.property(astra.FieldProperties.NAMEREF)) and field.type == astra.PropType.INT
            if add_to_required:
                required.append(field.name)
                        
        table_data["properties"] = fields_data
        table_data["required"] = required
        json_schema["properties"][table.name] = {"type": "array", "items": table_data}



    python_code = generate(
        json.dumps(json_schema),
        input_file_type=InputFileType.JsonSchema,
        output_model_type="pydantic_v2.BaseModel",
        target_python_version=PythonVersion.PY_310,
        formatters=[Formatter.BLACK, Formatter.ISORT],
        allow_population_by_field_name=True 

    )
    with open("models.py", "w", encoding="utf-8") as f:

        f.write(python_code)

    return {}

generate_model()
    
